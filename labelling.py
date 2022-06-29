import numpy as np

delta = 0.5

def get_vienna(filename):
    with open(filename) as f:
        return f.read()

def label_hairpins(vienna):
    l = len(vienna)
    output = np.zeros(l)
    r = 0.0

    for i, c in enumerate(vienna):
        if c == '(':
            r += delta
        elif c == ')':
            r -= delta

        output[i] = r

    return output

def save_hairpin_labeling(r_array):
    df = pd.DataFrame(r_array).T
    df.insert(0, "read-id", None)
    df.to_csv("hairpin_labeling.csv")
