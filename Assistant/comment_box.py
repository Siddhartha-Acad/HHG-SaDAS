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

def center_comment(text: str, width: int = 50, border: str = "~") -> str:
    line = f"# {border * (width - 2)}"
    middle = f"# {text.center(width - 4)} |"
    return "\n".join([line, middle, line])

def pretty_title(text: str, width: int = 50, fill: str = "~") -> str:
    return (f"{text}".center(width, fill)).replace(text, f": {text} :")


# comment_head_string = "Printing info"
# print(center_comment(comment_head_string))


comment_title_string = "HHG-SaDAS by Siddhartha Mithiya 2025"
print(pretty_title(comment_title_string))
