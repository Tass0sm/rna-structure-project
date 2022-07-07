def series_jerk(s):
    return s.diff().diff().fillna(value=0.0)

def series_jerk_points(s):
    jerk = series_jerk(s)

    return list(map(int, jerk[jerk != 0.0].index))
