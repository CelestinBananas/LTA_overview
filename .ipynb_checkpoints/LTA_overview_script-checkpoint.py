import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import itertools
from matplotlib.lines import Line2D

from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.visualization.wcsaxes import SphericalCircle



def split_saps(x):
    try:
        return x.split(",")
    except AttributeError:
        return None

def df_prep(datafile, to_keep, Debug = False):
    """
    Get relevant info for a selected list of projects as a new dataframe.

    Parameters
    ----------
    datafile : str
        The filename of the input database (i.e. with standardazed header structure)
    to_keep : list
        The list of projects to be selected for plotting
    Debug: bool, optional
        Print auxiliary info for verifications

    Returns
    -------
    result : dataframe
        Dataframe with the subset of info relevant for the seected projects
    """
    df_og = pd.read_csv(datafile, delimiter="\t", dtype ={'SAPS': str})
    df = df_og.copy()
#   df = df.drop(['SUCCESSOR_MOM_ID', 'SUCCESSOR_SAS_ID', 'SUCCESSOR_TYPE', 'STORAGE_MANAGER',
#                           'SPACE_USED', 'Size dysco compr [TB]',
#                           'Expected compression time if LBA per subband',
#                           'Expected compression time if HBA per subband',
#                           'Expected compression time per dataset' ], axis = 1)
    if Debug:
        print(f"Length of original df: {len(df.index)}")
    df = df[df["PROJECT"].isin(to_keep)]

    if Debug:
        print(f"Length using only projects in column G: {len(df.index)}")

    # A special case is represented by the project WTG-verification-DMO, for which main info is not available !
    if 'WTG-verification-DMO' in to_keep:
        print("WARNING project 'WTG-verification-DMO' is unavailable for plotting and has been automatically discarded from plots and dataframes.")
    df = df[df["PROJECT"] != 'WTG-verification-DMO']

    df = df.replace({'DEC': '0'}, '0,0')

    df['SAPS'] = df.SAPS.apply(lambda x: split_saps(x))

    #remove trailing commas
    df["RA"] = df.apply(lambda x: x.RA if x.RA[-1] != ',' else x.RA[:-1], axis = 1)
    df["DEC"] = df.apply(lambda x: x.DEC if x.DEC[-1] != ',' else x.DEC[:-1], axis = 1)

    #split observations with mulitple SAPS into list. Also need to deal with delimiter inconsistency ("," when alone, "." when in list of observations)
    df["RA"] = df.apply(lambda x: x.RA.split(',') if x.SAPS != ["0"] else [x.RA.replace(",", ".")], axis = 1)
    df["DEC"] = df.apply(lambda x: x.DEC.split(',') if x.SAPS != ["0"] else [x.DEC.replace(",", ".")], axis = 1)

    #remove weird extra 0s at the beginning of the RAs
    df["RA"] = df.apply(lambda x: x.RA if len(x.RA)==len(x.DEC) else x.RA[1:], axis = 1)

    #correct the SAPS that have missmatched numbers of SAPS compared to coords
    # TO BE DONE: VERIFY THIS ASSUMPTION
    df["SAPS"] = df.apply(lambda x: len(x.RA)*[np.nan] if x.SAPS == None else [i for i in range(len(x.RA))], axis = 1)


    # only exploded ra and dec as that is what is needed now but will probably have to come back
    # to this later if new info is required
    df = df.explode(["RA", "DEC", "SAPS"])
    df = df.reset_index(drop = True)

    return df


def overview_image(df, outname, plotProject = 'all', plotAntenna = 'both', legend_ncols = 2, Debug = False, testMode = False, plotAteam = False, exp_range = [- np.inf, np.inf]):
    """
    Plotting routine for a selected list of projects.

    Parameters
    ----------
    datafile : str
      The filename of the input database (i.e. with standardazed header structure)
    outname : str
      The prefix name for the generated outputs, e.g. plots
    plotProject : str
      The filter to select a subset of projects
    plotAntenna : str
      The filter for selecting and plottinf data of a specific antenna field
    legend_ncols : int
      The number of columns for the layout of the legend box
    Debug: bool, optional
      Print auxiliary info for verifications
    testMode : bool, optional
      Only limit the plotting to the first 100 datapoints
    plotAteam: bool, optional
        Plot the A team source circles with names, default False
    exp_range: list of ints, optional 
        first elemet is the minimum exposure in seconds and second is maximum exposure, default all exposures

    Returns
    -------
    result : plots
      Return specified plot
    """


    # cleanup Antenna column
    df["ANTENNA"] = df.apply(lambda x: x.ANTENNA_SET[:3], axis =1)

    #filter out duplicates in project, ra, dec and antenna
    df = df.drop_duplicates(subset = ["PROJECT", "RA", "DEC", "ANTENNA"])
    df = df[df["MIN_EXPOSURE"] < exp_range[1]]
    df = df[df["MIN_EXPOSURE"] > exp_range[0]]
    
    if Debug:
      print(f'Length without duplicates in Proj, Coords, Antenna:{len(df.index)}')

    #warn user if WTG selected as this was filtered out in df_prep
    if 'WTG-verification-DMO' in plotProject:
      print("WARNING project 'WTG-verification-DMO' is unavailable for plotting and has been automatically discarded from plots and dataframes.")

    #select only the projects that the user wants to plot creating new df "df_red"
    if plotProject == 'all':
      df_red = df
    else:
      #for 'all' keyword
      if isinstance(plotProject, str) and plotProject == 'all':
          plotProject = [plotProject]
          df_red = df[df["PROJECT"].isin(plotProject)]

      elif isinstance(plotProject, list):
          df_red = df[df["PROJECT"].isin(plotProject)]
      else:
          print("ERROR with plotProject input. Please use list of projects as strings or the 'all' keyword")

      if len(plotProject) != len(df_red["PROJECT"].unique()):
          print("!!! WARNING: Project list lengths differ, some projects may not be plotted. Please double check name and format. !!!")
          print(f"Length plotProject: {len(plotProject)}")
          print(f"Length of Projects to plot: {len(df_red['PROJECT'].unique())}")

    if Debug:
      print(f'Projects to plot: {df_red["PROJECT"].unique()}')

    #image stuff
    image = 'lambda_mollweide_haslam408_nofilt.fits'
    #     image = 'detthresh_P8R3_source_8years_PL22.fits'
    hdu = fits.open(image)[1]
    wcs = WCS(hdu.header)

    fig = plt.figure(figsize = (12, 5))
    fig.set_facecolor('white')

    ax = fig.add_subplot(projection = wcs, frame_class=EllipticalFrame)
    ax.imshow(hdu.data, vmin = 9e3, vmax = 1.5e5, cmap = 'hot')

    #hard code the A team sources in: [name, coords]
    A_team = [["Cas_A", "23h23m24s +58d48m54s"],
            ["Cyg_A", "19h59m28.3564s +40d44m02.096s"],
            ["Her_A", "16h51m08.024s +04d59m34.91s"],
            ["Tau_A", "05h34m31.94s +22d00m52.2s"],
            ["Ver_A", "12h30m49.4233s +12d23m28.043s"]]
    if plotAteam:
        for A_source in A_team:
            c = SkyCoord(A_source[1])
            pix = wcs.world_to_pixel(c.galactic)
#             ax.scatter(pix[0], pix[1], marker = '^', c = 'b', s = 20)
#             print(c.lon)
            sph = SphericalCircle((c.ra,c.dec), 4*u.deg, edgecolor = 'b', facecolor='none', linewidth = 1.5, transform= ax.get_transform('fk5'))
            ax.add_patch(sph)
            if A_source[0] == "Tau_A":
                ax.text(pix[0] - 450, pix[1] + 1, A_source[0], c = 'w')
            else:
                ax.text(pix[0] + 50, pix[1] + 1, A_source[0], c = 'w')

    #testmode to use a reduced number of points
    if testMode:
      df_to_run = df_red.iloc[:100]
    else:
      df_to_run = df_red

    #plot the points with appropirate transform and append proj name to the atenna lists. Later will use to
    # specify what projects are what antennas and count amount of points plotted
    HBA_proj = []
    LBA_proj = []

    #     markers = itertools.cycle(Line2D.markers)
    colours = itertools.cycle(get_cmap('tab20').colors)
    projs = df_to_run["PROJECT"].unique()

    def onpick3(event):
        ind = event.ind
        label = int(event.artist.get_label())
        print('Point selected:')
        print(f'\t Project code: {df_to_run["PROJECT"].loc[label]}')
        print(f'\t SAS ID: {df_to_run["SAS_ID"].loc[label]}')
        print(f'\t Antenna set: {df_to_run["ANTENNA_SET"].loc[label]}')
        print(f'\t Observation type: {df_to_run["OBSERVATION_TYPE"].loc[label]}')
        print(f'\t Min exposure: {df_to_run["MIN_EXPOSURE"].loc[label]} s')
        
    df_to_run["marker"] = df_to_run.apply(lambda x: "D" if x.ANTENNA == 'HBA' else "o", axis = 1)
    df_to_run["coords"] = df_to_run.apply(lambda x: SkyCoord(float(x.RA), float(x.DEC), unit = u.deg, frame='icrs'), axis = 1)
    df_to_run["pix"] = df_to_run.apply(lambda x: wcs.world_to_pixel(x.coords.galactic), axis = 1)
    
#     print(df_to_run["marker"])
#     pix0, pix1 = zip(*df_to_run["pix"])
#     ax.scatter(pix0[0],pix1[0], s = 20, label =df_to_run["PROJECT"], picker=True)
#     ax.scatter(pix0[1],pix1[1], s = 20, label =df_to_run["PROJECT"], picker=True)
    proj_colours = []
    for proj in projs:
      colour = next(colours)
      proj_colours.append([proj, colour])
      for index, row in df_to_run[df_to_run["PROJECT"] == proj].iterrows():
#           #make sure ra and dec in correct format for SkyCoord
# #           ra = float(row["RA"]) * u.deg
# #           try:
# #               dec = float(row["DEC"]) * u.deg
# #           except ValueError:
# #               print(index)
# #           coords = SkyCoord(ra, dec, unit = u.deg, frame='icrs')

#           #Plot only required antenna + append to list to know which antenna was used for label later (temporary: could maybe use df filtering?)
          if (row["ANTENNA"] == 'HBA') and ((plotAntenna == 'HBA') or (plotAntenna == 'both')):
#     #                 edgecolor = 'b'
#               sph_r = 0.5 * 3*u.deg
#               marker = "D"
              HBA_proj.append(row["PROJECT"])
#               HBA_proj.append(str(index))

#               pix = wcs.world_to_pixel(coords.galactic)
              ax.scatter(row["pix"][0], row["pix"][1], marker = row["marker"], color = colour, s = 20, label = index, picker = True)

          elif (row["ANTENNA"] == 'LBA') and ((plotAntenna == 'LBA') or (plotAntenna == 'both')):
#     #                 edgecolor = 'c'
#               sph_r = 0.5 * 13*u.deg
#               marker = "o"
              LBA_proj.append(row["PROJECT"])
#               LBA_proj.append(str(index))

#               pix = wcs.world_to_pixel(coords.galactic)
              ax.scatter(row["pix"][0], row["pix"][1], marker = row["marker"], color = colour, s = 20, label = index, picker = True)

    # shrink box by 20% to fit the legend outside the plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    #get and plot only the unique labels in the legend
    handles, labels = ax.get_legend_handles_labels()
#     print(handles)
#     print(df_to_run["PROJECT"].loc[0])
    labels = df_to_run["PROJECT"].loc[map(int, labels)]
#     unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    unique = [l for i, l in enumerate(labels) if l not in labels[:i]]
    unique.sort(key = lambda x:x[1])
#     print(unique)
#     a, b = zip(*unique)
#     print(b)
    

    total_numpoints = 0 #counter for number of plotted points
    #use a new list (unique2) to modify the legend labels
#     unique2 = []
#     for i in unique:
#       lst = list(i)

#       numpoints = HBA_proj.count(i[1]) + LBA_proj.count(i[1])
#       total_numpoints += numpoints

#       if (i[1] in HBA_proj) and (i[1] in LBA_proj):
#           lst[1] += f' #{numpoints} (HBA & LBA)'
#       elif i[1] in HBA_proj:
#           lst[1] += f' #{numpoints} (HBA)'
#       elif i[1] in LBA_proj:
#           lst[1] += f' #{numpoints} (LBA)'
#       else:
#           print(f'Errrrooorr: {i[1]}')

#       unique2.append(tuple(lst))
    unique2 = []
    for i in projs:
#       lst = list(i)

      numpoints = HBA_proj.count(i) + LBA_proj.count(i)
      total_numpoints += numpoints

      if (i in HBA_proj) and (i in LBA_proj):
          lst = f'{i} #{numpoints} (HBA & LBA)'
          marker = 'D'
      elif i in HBA_proj:
          lst = f'{i} #{numpoints} (HBA)'
          marker = 'D'
      elif i in LBA_proj:
          lst = f'{i} #{numpoints} (LBA)'
          marker = 'o'
      else:
          print(f'Errrrooorr: {i}')
       
      for j in proj_colours:
          if j[0] == i:
             c = j[1]
      handle = Line2D([0], [0], marker = marker, color = c,markerfacecolor = c, label = lst)
#       unique2.append(tuple(lst))
      unique2.append(handle)
    
    if Debug:
      print(f"Total points plotted: {total_numpoints}")

    # Put a legend to the right of the current axis
#     ax.legend(*zip(*unique2), loc='center left', bbox_to_anchor=(1.1, 0.5), ncol = legend_ncols)
    ax.legend(handles = unique2, loc='center left', bbox_to_anchor=(1.1, 0.5), ncol = legend_ncols)

    fig.suptitle(f"{outname}", y = 0.9)
    plt.subplots_adjust(top = 0.8, right = 0.5)
    fig.savefig(f"{outname}.png", bbox_inches = 'tight')
    fig.canvas.mpl_connect('pick_event', onpick3)
    plt.show()



def pick_custom_hit():
    # picking with a custom hit test function
    # you can define custom pickers by setting picker to a callable
    # function.  The function has the signature
    #
    #  hit, props = func(artist, mouseevent)
    #
    # to determine the hit test.  if the mouse event is over the artist,
    # return hit=True and props is a dictionary of
    # properties you want added to the PickEvent attributes



    fig, ax = plt.subplots()
    ax.set_title('custom picker for line data')
    line, = ax.plot(rand(100), rand(100), 'o', picker=line_picker)
    fig.canvas.mpl_connect('pick_event', onpick2)

########### Script ######################

Debug = False
testMode = False

#cp from column G of summary sheet
to_keep_raw = 'WTG-verification-DMO','DDT001','DDT002','DDT1_002','DDT2_001','DDT3_003','DDT4_002','DDT4_003','LC0_003','LC0_009','LC0_012','LC0_015','LC0_017','LC0_020','LC0_024','LC0_028','LC0_032','LC0_039','LC0_041','LC1_036','LC1_038','LC2_006','LC2_036','LC3_004','LC3_012','LC4_010','LC4_016','LC6_018','LC8_033','LC9_002'
# to_keep_raw = ["LC0_012"]
#run scripts
df_raw = df_prep('data/raw_data.tsv', to_keep_raw, Debug = Debug)
overview_image(df_raw, 'LTA_raw_data',plotProject ='all',plotAntenna = 'both', Debug = Debug, testMode = testMode, plotAteam=True, exp_range = [-np.inf, 14000])
