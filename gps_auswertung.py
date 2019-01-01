# -*- coding: utf-8 -*-
#!/usr/bin/python
import os
import sys
import numpy as np
import scipy.signal
import time, datetime
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
#import matplotlib.colors as colors
sys.path.append("/home/ronny/Documents/python/")        # path to srtm
#import pexpect
import srtm                     # see https://pypi.python.org/pypi/SRTM.py and https://github.com/tkrajina/srtm.py
#import gpxpy                   # Doesn't work as hoped, see below

# Reads one or several gpx files and creates a graph for elevation and speed in respect to time and distance.
#  The gpx files need to be given as parameter when running the script
# Also draws the track in an SRTM (elevation) map.
# If a file with the places exists, then the places which are on the track will added to the graph
#   Needs to contain the following tab-separated values: City Name, Radius of the city, Longitude, Latitude
# If a html file exists, then information can be added to that file
#   Needs to contain the following line after the line the data will be added: <!--Daten zufuegen-->

if len(sys.argv)==1:
    print("Start programm with parameter 'gpx-File'")
    exit(1)
parneupfadzeit=600      # Nach soviel s annehmen, dass neuer Pfad beginnt
parneupfaddist=5        # Nach sovielfach der median distanz zwischen 2 Datenpunkten, dass neuer Pfad beginnt
paranzbeg=15    #Ueberpruefen so vieler Datenpunke
parhoehbeg=30   #Höhendifferenzen zum Beginn größer als dieser herausschmeißen, da nach Einschalten noch zu ungenau
sigma = 5               # Sigmaclipping for abs(elevation), local speed
paranzhoeh=5            # Anzahl der Hoehen, über welche für Berechnung der Gesamthöhe gemittelt wird
paranzdist=15           # Lokale Geschwindigkeit über wie viele Datenpunkte berechnen
parrangeschw=[2,98]    #Min und Max zu plottende Geschwindigkeit in % aller lokalen Geschwindigkeiten

pdf_viewer  = ''#'okular'       # If empty the pdfs wo't be displayed
places_file = 'gaw_orte.dat'        # If it doesn't exists, places won't be read
html_file   = 'index.html'      # If it doesn't exists, than won't be processed
upload_data_command = 'gaw_upload.dat'      # If it doesn't exists, than won't be processed. Contains the upload command with a wildcart {0}

gradbog=np.pi/180


def sigma_clip(xarr, yarr, p_orders, sigma_l, sigma_h, repeats = 1):
    """
    Performs a sigmaclipping after fitting a polynomial against data
    :param xarr: 1d array with the x-values of the data
    :param yarr: 1d array with the data
    :p_orders: orders of the polynomial
    :sigma_l: Data off by this sigma are rejected on the lower side of the fit
    :sigma_h: Data off by this sigma are rejected on the higher side of the fit
    :repeats: redo the fit with the (cleaned) data how many times? 
    :return goodvalues: 1d array with True/False. The values which are inside the limits are True
    :return p: parameters of the last polynomial fit
    """
    xarr, yarr = np.array(xarr), np.array(yarr)
    if xarr.shape[0] == 0 or yarr.shape[0] == 0:
        logger('Warn: empty array for sigma clipping')
        return np.array([]), np.repeat([0], p_orders)
    elif xarr.shape[0] <> yarr.shape[0]:
        logger('Warn: got different sized arrays for sigma clipping. This is a programming error and should not happen')
        return [], np.repeat([0], p_orders)
    goodvalues = (yarr*0 == 0)
    old_values = [0,0]
    for i in range(repeats):
        poly = np.polyfit(xarr[goodvalues], yarr[goodvalues], p_orders)
        stddiff = np.std(yarr[goodvalues] - np.polyval(poly, xarr[goodvalues]), ddof=p_orders+1)
        diff = (yarr - np.polyval(poly, xarr))/stddiff
        diff[np.isnan(diff)] = 2*sigma_h                           # remove nans by setting it to the ourside value
        goodvalues = ( (diff >= -sigma_l) & (diff <= sigma_h) )    #average should be 0
        if stddiff == old_values[0] and poly[0] == old_values[1]:
            break
        old_values = [stddiff,poly[0]]
    #plot_img_spec.plot_spectra(np.array([xarr,xarr]),np.array([yarr,np.polyval(poly, xarr)]),['data','fit'], ['savepaths'], True, [0.06,0.92,0.95,0.06, 1.0,1.01], 'line {0}'.format(i))
    return goodvalues, poly

elevation_data = srtm.get_data()

orte, orte_namen = [], []
if os.path.isfile(places_file):
    file=open(places_file,'r')
    for line in file:
        if line.find('ä')>-1 or line.find('Ä')>-1 or line.find('ö')>-1 or line.find('Ö')>-1 or line.find('ü')>-1 or line.find('Ü')>-1 or line.find('ß')>-1:
            print("Sonderzeichen:",line[:-1])
        line=line[:-1].split('\t')
        if len(line) < 4:
            print "Eintrag zu kurz:",line
            exit(1)
        orte.append([float(line[3]),float(line[2]),float(line[1])/40000*2*np.pi])    # Breite, Länge, Radius[rad]
        orte_namen.append(line[0])      # Name Ort
    file.close
    print "Orte eingelesen"
    
#print('reading gpx data')
orte = np.array(orte)
gpx_files=[]
daten=[]        # index of filename (pos in gpx_files)
neue_dat=[]     # remove
for gpxdat in sys.argv[1:]:
    ### schauen, ob Name richtig
    if gpxdat[-3:]==".gz":
        gpxdat=gpxdat[:-3]
    if gpxdat[-4:]==".gpx":
        gpxdat=gpxdat[:-4]
    if gpxdat[-1]==".":
        gpxdat=gpxdat[:-1]
    temp=os.popen('ls '+gpxdat+'.gpx*').readlines()
    if len(temp)==0:
        print "keine passende Datei:",gpxdat+'.gpx*'
        exit(1)
    if len(temp)>1:
        print "mehrere passende Dateien:",gpxdat+'.gpx*'
        exit(1)
    if temp[0].find('.gz')>0:
        os.system('gunzip '+temp[0])
    gpx_files.append(gpxdat)
    help=0
    neue_dat.append(len(daten))
    #gpx_data = gpxpy.parse(open(gpxdat+'.gpx','r'))            # Replace later? It doesn't work with my files
    file=open(gpxdat+'.gpx','r')
    for line in file:
        if line.find('<trkpt')==0:
            line=line.split('"')
            for k in [1,3]:                 # latitude, longitude
                line[k]=float(line[k])
            #index of filename,zeit in s,lat,lon,hoehe,dist,hoehenunterschied, geschw lokal, geschw ges, index naechster ort in orte/orte_name, Distanz zum naechsten Ort, Hoehe gefiltert, Elevation from SRTM
            new_data = [len(gpx_files)-1, np.nan, line[1], line[3], np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]   
            dist = np.arccos(np.sin(gradbog*line[1])*np.sin(gradbog*orte[:,0]) + np.cos(gradbog*line[1])*np.cos(gradbog*orte[:,0])* np.cos(gradbog*(line[3]-orte[:,1])) )
            min_dist = np.argmin(np.abs(dist))
            if dist[min_dist] < orte[min_dist,2]:
                new_data[9], new_data[10] = min_dist, dist[min_dist]
            help=1        #Wenn erster Trackpoint gefunden, dann alle lesen
            continue    #zur naechsten Zeile gehen, welche die Hoehe enthaelt
        if help==0:    #nicht in Trackpointumgebung
            continue
        if line.find('<ele>')==0:    #Hoehe hinzufüegen
            line=line.split('</ele>')
            new_data[4]=float(line[0][5:])
        elif line.find('<time>')==0:    #Zeit hinzufügen
            line = line.split('</time>')
            line = line[0].split('<time>')
            line = datetime.datetime.strptime(line[1], '%Y-%m-%dT%H:%M:%SZ')
            new_data[1] = time.mktime(line.timetuple())
        elif line.find('<type>')==0:        #nur Distanz-Positionen zulassen, wenn nicht vorhanden, dann egal
            line=line.split('<type>')
            line=line[1].split('</type>')
            if line[0].find('D')==-1 and line[0].find('V')==-1 and line[0].find('T')==-1:
                help = 0
        elif line.find('</trkpt>')==0:
            daten.append(new_data)
            help=0        #aktueller Trackpoint abgearbeitet.
    file.close()
    
daten = np.array(daten)
datens = daten.shape
print("File/Files read")

# Get the SRTM elevation
print('Get the SRTM elevations')
for i in range(datens[0]):
    #elevation_data.add_elevations(gpx_data)                    # How do I get the data out?
    #print elevation_data.get_elevation(50.8682, 7.1377)        # Does work, but only with floats, not lists
    daten[i,12] = elevation_data.get_elevation(daten[i,2], daten[i,3])

# Schlechte Hoehenmessung am Anfang aussortieren
new_file = np.where(daten[1:,0] - daten[:-1,0] <> 0)[0]
for j in new_file:    #fuer jede gpx-Datei
    badvalues = np.where(abs(daten[j:j+paranzbeg,4]-daten[j+paranzbeg,4]) > parhoehbeg)[0]     # Hoehenunterschied am Beginn Datei zu groß
    daten[badvalues+j,4:8] = np.nan         # Hoehe, Dist, Hoehenunterschied, Geschw

# Distance between consecutive points
dist = np.sin(gradbog*daten[1:,2])*np.sin(gradbog*daten[:-1,2]) + np.cos(gradbog*daten[1:,2])*np.cos(gradbog*daten[:-1,2])*np.cos(gradbog*(daten[1:,3]-daten[:-1,3]))
dist[dist > 1] = 1
dist[dist < -1] = -1
dist = np.arccos(dist) *6378*1000   # Erdradius in m
# Elevetion difference
hoehdiff = daten[1:,4]-daten[:-1,4]
med_dist = np.nanmedian(dist)
# Cleaning zu viel Zeit, zu viel Entfernung, neue gpx datei
new_dataset = np.where((daten[1:,1] - daten[:-1,1] > parneupfadzeit) | (dist > med_dist*parneupfaddist) | (daten[1:,0] - daten[:-1,0] <> 0) )[0]     # conditions conected by or
dist[new_dataset] = 0           #np.nan
hoehdiff[new_dataset] = np.nan
# Cleaning too big elevation steps
good_values, poly = sigma_clip(hoehdiff, np.abs(hoehdiff), 0, sigma, sigma, repeats=1)
hoehdiff[~good_values] = np.nan
good_values = np.insert(good_values, 0, False)
daten[~good_values,4] = np.nan
print('Excluded because of elevation: {0} of {1}'.format(len(daten[~good_values,4]), datens[0]))
# Hinzufuegen
daten[1:,5] = dist
daten[1:,6] = hoehdiff

# Medfilter for the elevations
new_dataset = np.append(new_dataset, datens[0] )
new_dataset = np.insert(new_dataset, 0, 0 )

for i in range(len(new_dataset)-1):
    mindat, maxdat = new_dataset[i], new_dataset[i+1]
    daten[mindat:maxdat,11] = scipy.signal.medfilt(daten[mindat:maxdat,4], paranzhoeh)  #
    daten[mindat:min(mindat+paranzhoeh,maxdat),11] = np.nan
    daten[max(maxdat-paranzhoeh,mindat):maxdat,11] = np.nan
    # Anstelle des Median Filters: Fit Polynom und dann  Sigmaclipping? But only 5 points?

med_dist = np.nanmedian(dist)
times = daten[1:,1] - daten[:-1,1]
times[times < 1] = 1            # At least 1s
med_geschw = [np.nanmedian(dist / times), np.nanstd(dist / times, ddof=1) ]
print('Median Werte zwischen aufeinanderfolgenden Datenpunkten: Entfernung: {0}m, Geschwindigkeit: {1}m/s (+-{2})'.format(round(med_dist,2), round(med_geschw[0],2), round(med_geschw[0],1) ))
pausendauer=med_dist/0.5    #Pause an Ampel oder aehnlich 0.5m/s==2km/h

temphoehe = daten[1:,11] - daten[:-1,11]            # Medianfiltered Hoehendifferenz zwischen 2 aufeinanderfolgenden Punkten, nan at data gaps
temphoehe = temphoehe[~np.isnan(temphoehe)]
hoehe = [np.sum(temphoehe[temphoehe > 0]), np.sum(temphoehe[temphoehe < 0])]    # Bergauf, Bergab
hoehen = [np.nanmin(daten[:,4]), np.nanmax(daten[:,4])]                         # Min Hoehe, Max Hoehe
hoehen.append((hoehen[1] + hoehen[0])/2)                                        # median elevation
hoehen.append((hoehen[1] - hoehen[0]))                                        # Elevation difference
slow = (dist/times > pausendauer) | (dist/times < med_geschw[0]/10)
times[slow] = dist[slow]/med_geschw[0]                                           # med time passed when pausing or too slow
times[np.isnan(dist)] = 0
dauer = np.nansum(times)

for i in range(paranzdist+1,datens[0]-paranzdist/2):
    daten[i,8] = np.nansum(dist[:i]) / np.nansum(times[:i]) * 3.6                                    # Global speed up to here
    daten[i,7] = np.nansum(dist[i-1-paranzdist/2:i+paranzdist/2]) / np.nansum(times[i-1-paranzdist/2:i+paranzdist/2]) * 3.6   # Local speed here
good_values, poly = sigma_clip(daten[:,7], daten[:,7], 0, sigma, sigma, repeats=1)
daten[~good_values,7] = np.nan
speeds = [np.nanpercentile(daten[:,7],parrangeschw[0]), np.nanpercentile(daten[:,7],parrangeschw[1])]             # min and max speed to plot

# Remove places which are doubled
delort=[]
with_places = np.where(~np.isnan(daten[:,10]))[0]
for i in with_places:
    for j in range(min(i-1,20)):                #Bestimmen des kleinsten Abstands zum Ort
        if daten[i,9] == daten[i-(j+1),9]:        # same place
            if daten[i,10] < daten[i-(j+1),10]:   # smaller distance
                daten[i-(j+1),10] = np.nan         # remove the previous place
                daten[i-(j+1),9] = np.nan          
            else:
                delort.append(i)                    # remove this place later
                break
daten[delort,10]=np.nan
daten[delort,9]=np.nan


# Testing:
"""
geo_elevation_data = srtm.get_data()
image = geo_elevation_data.get_image((500, 500), (45, 45.5), (13.5, 14), 300)
# the image s a standard PIL object, you can save or show it:
image.show()
image_arr = geo_elevation_data.get_image((500, 500), (45, 45.5), (13.5, 14), 300, mode='array')
elev_range = [np.nanmin(image_arr), np.nanmax(image_arr)]
print elev_range"""

# Read the data for the map, note that latitude is on x-axis, because python is plotting in that way
# SRTM1 (for US territories) with 1 arcsec resolution), STRM3 (for world with 3arcsec resolution), on both axis without cos(lat) correction. This means the longitude is with better resolution than lattitude
# using lats[0] for the correction of longitude as this is the zeropoint of the image
lats  = [np.nanmin(daten[:,2]), np.nanmax(daten[:,2])]
longs = [np.nanmin(daten[:,3]), np.nanmax(daten[:,3])]
dlats, dlongs = max(1./30., lats[1] - lats[0]), max(1./30., longs[1] - longs[0])
lats = [lats[0] - 0.1 * dlats, lats[1] + 0.1 * dlats]               # add some area around 
longs = [longs[0] - 0.1 * dlongs, longs[1] + 0.1 * dlongs]          # add some area around 
factor = 3 * np.cos(gradbog*lats[0])          # 1 for SRTM1 and 3 for SRTM3, corrected with cos(lat) to use all resolution, also in longitude.
#factor = 111000 / (max(med_dist,100)/2)                              # 1deg is about 111km, sampling over half the median distance/90m
x_arr = (np.round((daten[:,2] - lats[0]) * 3600 / factor) ).astype(int)                                  # Convert measured coordinates into px
y_arr = (np.round((daten[:,3] - longs[0]) * np.cos(gradbog*lats[0]) * 3600 / factor) ).astype(int)    # Convert measured coordinates into px
#y_arr = (np.round((daten[:,3] - longs[0]) * 3600 / factor) ).astype(int)    # Convert measured coordinates into px
size_x = int((lats[1] - lats[0]) * 3600 / factor) + 1
size_y = int((longs[1] - longs[0]) * np.cos(gradbog*lats[0]) * 3600 / factor) + 1
size_km = [111 * (lats[1] - lats[0]), 111 * (longs[1] - longs[0]) * np.cos(gradbog*lats[0])]
#size_y = int((longs[1] - longs[0]) * 3600 / factor) + 1
#image = elevation_data.get_image((size_y, size_x), (lats[0]-0.5, lats[1]+.5), (longs[0]-.5, longs[1]+.5), 300)   # (1000, 1000): image size; 300: maximale hoehe fuer hellstes gruen
# the image s a standard PIL object, you can save or show it:
#image_arr = elevation_data.get_image((size_y, size_x), (lats[0], lats[1]), (longs[0], longs[1]), 300, mode='array')
#image = elevation_data.get_image((size_y, size_x), (lats[0], lats[1]), (longs[0], longs[1]), 300)
#image.show()

# Fill the image array with data
image_arr = np.empty((size_x,size_y))       # y and x are flipped because of the way python is plotting
image_arr.fill(np.nan)
latitude = np.linspace(lats[0], lats[1], size_x)
longitude = np.linspace(longs[0], longs[1], size_y)
for i in range(len(latitude)):
    for j in range(len(longitude)):
        elevation = elevation_data.get_elevation(latitude[i], longitude[j])
        image_arr[i,j] = elevation
# Get the elevation borders
if len(image_arr[~np.isnan(image_arr)]) > 0:
    elev_range = [np.nanmin(image_arr), np.nanmax(image_arr)]
else:
    elev_range = [-1, -1]
    print('Problem: no SRTM data available')
delev_range = elev_range[1] - elev_range[0]
# Plot the path
daten[:,13] = image_arr[x_arr, y_arr]                   # get the data from the image
image_arr[x_arr, y_arr] = elev_range[1] + 0.1 * delev_range # replace the path with the brightest data


print "Finished calculations, started plotting..."

if elev_range <> [-1, -1]:
    # Plotting the image with the path
    fontscale = 2.5
    title = 'STRM (elevation) map of tour'
    for i in gpx_files:
        title += ' {0},'.format(i)
    fig, frame = plt.subplots(1, 1)
    fig.set_size_inches(32, 32)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    cframe = frame.imshow(image_arr, origin='lower', cmap='YlGn_r',vmin=elev_range[0], vmax=elev_range[1] + 0.1 * delev_range, interpolation='none', extent=[0,size_km[1],0,size_km[0]]) #cmap='gray'
    cframe.axes.tick_params(rotation=0, labelsize=12*fontscale)
    frame.set_xlabel('Distance [km], Longitude: {0} ... {1}'.format(round(longs[0],3), round(longs[1],3)), fontsize=14*fontscale)
    frame.set_ylabel('Distance [km], Latitude: {0} ... {1}'.format(round(lats[0],3), round(lats[1],3)), fontsize=14*fontscale)
    frame.set_title(title[:-1], fontsize=16*fontscale)
    cbar = plt.colorbar(cframe,fraction=0.024, pad=0.02)
    cbar.set_label('Elevation [m] ({0} m - {1} m)'.format(int(round(elev_range[0])), int(round(elev_range[1]))), fontsize=14*fontscale)
    plt.locator_params(nbins=20)            # More Tick labels at the axes
    plt.savefig(gpx_files[0]+'_map.png', bbox_inches='tight')
    plt.close()
#exit(1)

# Plotting the profile
daten[:,5] = np.nancumsum(daten[:,5])/1000               # Distance to this point
distance = [np.nanmin(daten[:,5]), np.nanmax(daten[:,5]), np.nanmax(daten[:,5])-np.nanmin(daten[:,5]) ] # min, max, dx
timedelta = daten[-1,1] - daten[0,1]
if timedelta > 86460:                   # More than one day
    timestr = "%d days, %H:%M"
else:
    timestr = "%H:%M"
legende = 'Total: Elevation up: {0} m, Elevation down: {1} m,\nDistance: {2} km,\nTime: {3}, Time moving: {4}'
legende = legende.format(int(hoehe[0]), int(hoehe[1]), round(distance[1],2),
                        time.strftime(timestr, time.gmtime(timedelta)), time.strftime("%H:%M", time.gmtime(np.nansum(times))) )
fig, frame1 = plt.subplots(1, 1)
fig.set_size_inches(16.2, 10)
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
frame2 = frame1.twinx()
frame1.plot(daten[:,5], daten[:,12], label='Elevation SRTM', linestyle='-',marker='', color='brown')
frame1.plot(daten[:,5], daten[:,4], label='Elevation GPS', linestyle='',marker='+', color='r')
frame1.plot(daten[:,5], daten[:,11], label='Elevation GPS smoothed ({0} points)'.format(paranzhoeh), linestyle='',marker='.', color='purple', markersize=1)
#frame1.plot(daten[:,5], daten[:,13], label='Elevation SRTM (image)', linestyle='-',marker='', color='orange')      # was only added for testing
frame2.plot(daten[:,5], daten[:,7], label='local speed ({0} points)'.format(paranzdist), linestyle='',marker='x', color='chartreuse')
frame2.plot(daten[:,5], daten[:,8], label='speed up to here', linestyle='',marker='*', color='dodgerblue')
with_places = np.where(~np.isnan(daten[:,10]))[0]
lastplace = -1E6
ortschaften = ''
for i in with_places:
    if daten[i,5] - lastplace < distance[2]*0.02:
        continue
    lastplace = daten[i,5]
    ortschaften += ' - ' + orte_namen[int(daten[i,9])]
    for mult in range(1,1000):                                  # Increase the search area in case no elevation is available in the smallest search area
        daten_y = daten[max(0,i-int(0.008*datens[0]*mult)):min(datens[0],i+int(0.008*datens[0]*mult)+1),4]              # Elevation around the point
        if len(daten_y[~np.isnan(daten_y)]) >= 3:
            break
    # print daten[i,11], daten[i,4], daten[i,12], daten[i,13], hoehen, np.nanpercentile(daten_y, 10)          # Check why positions of places are so far off
    if np.nanmedian([daten[i,11], daten[i,4], daten[i,12], daten[i,13]]) > hoehen[2]:         # Bigger than middle height, use different elevations in case one is NaN (alternatively, use area [i-x:i+x+1]?)
        verticalalignment = 'top'
        y_pos = np.nanpercentile(daten_y, 10) - 0.05 * hoehen[3]                              # plot at bottom
    else:
        verticalalignment = 'bottom'                       # plot at top
        y_pos = np.nanpercentile(daten_y, 90) + 0.05 * hoehen[3]
    frame1.text(daten[i,5], y_pos, orte_namen[int(daten[i,9])],
               horizontalalignment='center', verticalalignment=verticalalignment,
               rotation=90, color='k', zorder=15, fontsize=15)
for i in range(0, datens[0], max(1,int(0.02*datens[0])) ):                    # Plot here to make sure the times are at the bottom of the scale
    frame2.text(daten[i,5], speeds[0], datetime.datetime.fromtimestamp(daten[i,1]).strftime("%H:%M"),
               horizontalalignment='center', verticalalignment='bottom',
               rotation=90, color='k', zorder=20, fontsize=10)
frame1.set_zorder(frame2.get_zorder()+1)                    # put frame1 in front of frame2
frame1.patch.set_visible(False)                             # hide the canvas
frame1.set_xlabel('Distance [km]', fontsize=15)
frame1.set_ylabel('Elevation [m]', fontsize=15)
frame2.set_ylabel('Speed [km/s]', fontsize=15)
frame1.set_title(legende, fontsize=17)
frame1.legend(loc='upper left',  bbox_to_anchor=(-0.05, 1.1), fontsize=15)
frame2.legend(loc='upper right', bbox_to_anchor=(1.05, 1.1), fontsize=15)
plt.axis([distance[0]-distance[2]*0.01, distance[1]+distance[2]*.01, speeds[0], speeds[1]])
plt.locator_params(nbins=20)
plt.savefig(gpx_files[0]+'.pdf', bbox_inches='tight')
plt.savefig(gpx_files[0]+'.png', bbox_inches='tight')
plt.close()
#print("\nMehrere PDFs zusammenfuegen: pdftk tour*.pdf cat output ziel.pdf\n")
if pdf_viewer <> '':
    os.system('{1} {0}.pdf &'.format(gpx_files[0], pdf_viewer))


############# Einbinden in html-Datei
subfolder = ''              # can be a folder or a complete address, e.g. http://h2324143.stratoserver.net/~ronnyabroad/gps/
if os.path.isfile(html_file):
    html=[]
    adden=1
    file=open(html_file,'r')
    for line in file:
        html.append(line)
        if line.find(gpx_files[0])>0:       # only add the file, if it's not already in the html table
            adden=0
    file.close()
    if adden==1:
        href=''
        for i in gpx_files:
            href=href+' <a href="{1}{0}.gpx">gpx</a>'.format(i, subfolder)
        file=open('index.html','w')
        for line in html:
            if line.find('<!--Daten zufuegen-->')>-1:
                file.write('  <tr>\n')
                # name:
                file.write('    <td width="68%" align="left"><font size=-1><b>'+gpx_files[0]+'</b>'+ortschaften+'</font></td>\n')
                # pdf, png, map.png
                file.write('    <td width="6%" align="center"><font size=-1> <a href="{1}{0}.pdf">pdf</a>, <a href="{1}{0}.png">png</a</font></td>\n'.format(gpx_files[0], subfolder))
                if elev_range <> [-1, -1]:
                    file.write('    <td width="6%" align="center"><font size=-1> <a href="{1}{0}_map.png">map-png</a</font></td>\n'.format(gpx_files[0], subfolder))
                else:                                               # if no SRTM data is available
                    file.write('    <td width="6%" align="center"><font size=-1> <a href="{1}{0}_map.png"></a</font></td>\n'.format(gpx_files[0], subfolder))
                # html: 
                temp=href.replace('gpx.gz','html')
                temp=temp.replace('gpx','html')
                file.write('    <td width="8%" align="center"><font size=-1>'+temp+'</font></td>\n')
                # gpx:
                file.write('    <td width="6%" align="center"><font size=-1>'+href+'</font></td>\n')
                # kml:
                file.write('    <td width="6%" align="center"><font size=-1>'+href.replace('gpx','kml')+'</font></td>\n')            
                file.write('  </tr>\n')
            file.write(line)
        file.close()
#for i in gpx_files:
    #os.system('sed -i \'s/amp" type=/amp;key=ABQIAAAA73cUYNc6rT4y9mjKEclHShSBM92zzWw9S0fp7QW1vjI8ZfmM0BS9Ga50GiYD85Iw7G50vAEIH2ow_w" type=/g\' '+ i+'.html')
    #print 'sed -i \'s/amp" type=/amp;key=ABQIAAAA73cUYNc6rT4y9mjKEclHShSBM92zzWw9S0fp7QW1vjI8ZfmM0BS9Ga50GiYD85Iw7G50vAEIH2ow_w" type=/g\' '+ i+'.html'

############# Hochladen auf Server

if os.path.isfile(upload_data_command):
    text='index.html Touren.kmz ../wichtigeProtokolle/fahrrad.ods {0}.pdf {0}.png {0}_map.png'.format(gpx_files[0])
    for j in gpx_files:
        text=text+' '+j+'.html '+j+'.gpx '+j+'.kml'
    
    scp_commands = open(upload_data_command,'r').readlines()
    for scp_command in scp_commands:
        if len(scp_command) < 5:
            continue
        print('Copying data')
        scp_command = scp_command[:-1].format(text)      # without \n but with replaced text
        #print scp_command
        scp = os.system(scp_command)
        if scp > 0:
            print("Probleme mit scp, Exitstatus: {0}\n\nIf necessary rerun\n {1}".format(scp, scp_command))
        else:
            print('Finished this data transfer')






