
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

