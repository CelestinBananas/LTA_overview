# LTA_overview

Version 0 :
- [x] Split HBA LBA, project, debug mode: include verification checks
- [x] TBD: add A-team sources location on the skymap view; add warning for filtering out the WTG-verification-DMO
- [x] Add iteration over markers or colours for each prjoect

Version 1 (Tuesday 25/01/22):
- [x] Include exposures -filter in plot function for exposure range, default all 
- [x] Circle for A-team and make optional
- [x] Clickable markers - proj, SASID, antenna, ant filter
- [x] Make background optional
- [ ] Implement interactive plots on notebook
- [ ] Add Fermi gamma ray,  2MASS and TGSS maps https://irsa.ipac.caltech.edu/Missions/2mass.html https://vo.astron.nl/tgssadr/q_fits/imgs/form


Version 2:
- [ ] 90% Time and smearing
- [ ] lb data sets (output) (low prio)
---
Other Notes:
- Reevaluate the usefulness of the to_keep list, is this not similar/or made redundant with the plotprojects specification of the plot_overview function ?
- Filtering unique ra and dec and antenna is temp, later will take into account 
- Will eventually run with Julich data (will require extra input param) and add if project/data is in both datbases

Notes from 12/01:
- clickable markers
- interactive!
