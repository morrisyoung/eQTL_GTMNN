import numpy as np




def exact_mc_perm_test(xs, ys, nmc):
    n, k = len(xs), 0
    diff = np.abs(np.mean(xs) - np.mean(ys))
    zs = np.concatenate([xs, ys])
    print zs

    for j in range(nmc):
        np.random.shuffle(zs)
        print zs
        k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
    return float(k) / nmc






if __name__ == "__main__":


	xs = np.array([12.6, 11.4, 13.2, 11.2, 9.4, 12.0])
	ys = np.array([16.4, 14.1, 13.4, 15.4, 14.0, 11.3])

	print exact_mc_perm_test(xs, ys, 30000)



