"""
File: manim_colors.py
Project: HHG-SaDAS

Code Description:
    | This file defines a collection of color palettes for plotting
    | and visualization within the HHG-SaDAS project.
    |
    | These palettes provide a consistent, visually appealing color scheme
    | for figures and data visualization, including those used in the thesis.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
-
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


import numpy as np

#          ~: A :~    ~: B :~    ~: C :~    ~: D :~    ~: E :~
blue   = ['#C7E9F1', '#9CDCEB', '#58C4DD', '#29ABCA', '#1C758A']     # ~: 0 :~
green  = ['#C9E2AE', '#A6CF8C', '#83C167', '#77B05D', '#699C52']     # ~: 1 :~
yellow = ['#FFF1B6', '#FFEA94', '#FFFF00', '#F4D345', '#E8C11C']     # ~: 2 :~
gold   = ['#F7C797', '#F9B775', '#F0AC5F', '#E1A158', '#C78D46']     # ~: 3 :~
red    = ['#F7A1A3', '#FF8080', '#FC6255', '#E65A4C', '#CF5044']     # ~: 4 :~
maroon = ['#ECABC1', '#EC92AB', '#C55F73', '#A24D61', '#94424F']     # ~: 5 :~
purple = ['#CAA3E8', '#B189C6', '#9A72AC', '#715582', '#644172']     # ~: 6 :~

des_col_1 = ['m', '#398C46', '#2D70B4', '#FA7E1A', '#C84440', '#6042A6']
des_col_2 = ['mediumslateblue', '#398C46', '#2D70B4', '#FA7E1A', '#C84440', '#6042A6']
named_color_1 = ['silver', 'crimson', 'mediumspringgreen', '#C71585', 'aquamarine', 'green', '#A87945', '#FF9CA0', 'orangered']
named_color_2 = ['silver', 'crimson']

A = ['#C7E9F1', '#CAA3E8', '#C9E2AE', '#F7C797', '#FFF1B6', '#F7A1A3']
B = ['#9CDCEB', '#B189C6', '#A6CF8C', '#F9B775', '#FFEA94', '#FF8080']
C = ['#58C4DD', '#9A72AC', '#83C167', '#F0AC5F', '#FFFF00', '#FC6255']
D = ['#29ABCA', '#715582', '#77B05D', '#E1A158', '#F4D345', '#E65A4C']
E = ['#1C758A', '#644172', '#699C52', '#C78D46', '#E8C11C', '#CF5044']

teal_c = '#5CD0B3'
sweet_green = '#83C167'
git_green = '#70BF41'
HTML_sweet_green = '#16a085'


C_L = ['#58C4DD', '#9A72AC', '#83C167', '#F0AC5F', '#E8C11C', '#FC6255']
C_L_reordered = ['#58C4DD', '#9A72AC', '#83C167', '#F0AC5F', '#3A7AB7']     #3A7AB7 <- darker version of the last one (6A9FD1).
GnuPlot_colors = ['indigo', '#009e73', '#56b4e9', '#e69f00', '#0072b2', '#e51e10', '#f0e442']


# ~~~~: https://github.com/karthik/wesanderson/tree/master :~~~~~
wes_palettes = {
    "BottleRocket1": np.array(["#A42820", "#5F5647", "#9B110E", "#3F5151", "#4E2A1E", "#550307", "#0C1707"]),
    "BottleRocket2": np.array(["#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"]),
    "Rushmore1": np.array(["#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F"]),
    "Rushmore": np.array(["#E1BD6D", "#EABE94", "#0B775E", "#35274A", "#F2300F"]),
    "Royal1": np.array(["#899DA4", "#C93312", "#FAEFD1", "#DC863B"]),
    "Royal2": np.array(["#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089"]),
    "Zissou1": np.array(["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]),
    "Zissou1Continuous": np.array(["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05", "#EF5703", "#F11B00"]),
    "Darjeeling1": np.array(["#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"]),
    "Darjeeling2": np.array(["#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000"]),
    "Chevalier1": np.array(["#446455", "#FDD262", "#D3DDDC", "#C7B19C"]),
    "FantasticFox1": np.array(["#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"]),
    "Moonrise1": np.array(["#F3DF6C", "#CEAB07", "#D5D5D3", "#24281A"]),
    "Moonrise2": np.array(["#798E87", "#C27D38", "#CCC591", "#29211F"]),
    "Moonrise3": np.array(["#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"]),
    "Cavalcanti1": np.array(["#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15"]),
    "GrandBudapest1": np.array(["#F1BB7B", "#FD6467", "#5B1A18", "#D67236"]),
    "GrandBudapest2": np.array(["#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4"]),
    "IsleofDogs1": np.array(["#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680"]),
    "IsleofDogs2": np.array(["#EAD3BF", "#AA9486", "#B6854D", "#39312F", "#1C1718"]),
    "FrenchDispatch": np.array(["#90D4CC", "#BD3027", "#B0AFA2", "#7FC0C6", "#9D9C85"]),
    "AsteroidCity1": np.array(["#0A9F9D", "#CEB175", "#E54E21", "#6C8645", "#C18748"]),
    "AsteroidCity2": np.array(["#C52E19", "#AC9765", "#54D8B1", "#b67c3b", "#175149", "#AF4E24"]),
    "AsteroidCity3": np.array(["#FBA72A", "#D3D4D8", "#CB7A5C", "#5785C1"])
}
