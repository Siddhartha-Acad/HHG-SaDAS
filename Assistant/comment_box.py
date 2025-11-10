"""
File: comment_box.py
Project: HHG-SaDAS
Code Description:
    | Utility script to generate visually centered comment boxes.

Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""

# for function_bank.py : comment box size = 60
# for normal codes     : comment box size = 50

comment_symbol = ''            # comment symbol: '#', '!', '\\'
def center_comment(text: str, width: int = 70, border: str = "~") -> str:
    line = f"{comment_symbol} {border * (width - 2)}"
    middle = f"{comment_symbol} {text.center(width - 4)} |"
    return "\n".join([line, middle, line])

def pretty_title(text: str, width: int = 40, fill: str = "~") -> str:
    line = f": {text} :"
    return f"{comment_symbol} " + line.center(width - 1, fill)


# comment_head_string = "HHG_spectra.py"
# print(center_comment(comment_head_string))

comment_title_string = "LASER info"
print(pretty_title(comment_title_string))
