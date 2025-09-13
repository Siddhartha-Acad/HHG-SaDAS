"""
File: Time_conversion.py
Project: HHG-SaDAS

Code Description:
    | This module provides utility functions to convert time in seconds
    | into minutes, hours, or days with appropriate formatting.
    |
    | Functions included:
    |   • secs_to_min_sec(seconds): Returns (minutes, seconds)
    |   • seconds_to_hours_minutes_seconds(seconds): Returns (hours, minutes, seconds)
    |   • secs_to_hr_min_sec(seconds): Returns (hours, minutes, seconds)
    |   • secs_to_day_hr_min(total_seconds): Returns (days, hours, minutes)
    |
    | These functions are useful for reporting elapsed time or converting
    | computational durations into readable formats.


Author: Siddhartha Mithiya
Affiliation: Indian Institute of Technology (IIT) Mandi
License: MIT License
Repository: https://github.com/Siddhartha-Acad/HHG-SaDAS.git

--------------------------------------------------------------------------------
Notes:
- Will be used as a module to convert measured execution time of programmes.
- This file is part of the HHG-SaDAS package, developed during my MS(R) thesis:
  "Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60".
--------------------------------------------------------------------------------
"""


def secs_to_min_sec(seconds):
    minutes = seconds // 60
    seconds %= 60
    return int(minutes), int(seconds)
def seconds_to_hours_minutes_seconds(seconds):
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return int(hours), int(minutes), int(seconds)
def secs_to_hr_min_sec(seconds):
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return int(hours), int(minutes), int(seconds)
def secs_to_day_hr_min(total_seconds):
    days = total_seconds // 86400
    hours = (total_seconds % 86400) // 3600
    minutes = (total_seconds % 3600) // 60
    return int(days), int(hours), int(minutes)

