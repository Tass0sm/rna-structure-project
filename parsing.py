def parse_until_zero(s, i):
    while not s.iloc[i] == 0.0:
        i += 1

        if i >= len(s) - 1:
            break

    return i

def parse_backward_until_zero(s, i):
    while not s.iloc[i] == 0.0:
        i -= 1

        if i <= 0:
            break

    return i

def parse_until_positive(s, i):
    while not s.iloc[i] > 0:
        i += 1

        if i >= len(s) - 1:
            break

    return i

def parse_until_negative(s, i):
    while not s.iloc[i] < 0:
        i += 1

        if i >= len(s) - 1:
            break

    return i

def parse_backward_until_positive(s, i):
    while not s.iloc[i] > 0:
        i -= 1

        if i <= 0:
            break

    return i + 1

def parse_loop(s, i):
    end = parse_until_negative(s, i)
    start = parse_backward_until_positive(s, end)

    i = end

    return i, (start, end)

def parse_top_loop(s, i):
    start = parse_until_positive(s.diff(), i)
    end = parse_until_zero(s, start)

    loop_start = int(s.iloc[start:end].idxmax()) + 1
    loop_end = parse_until_negative(s.diff(), loop_start)

    return end, (loop_start, loop_end)
