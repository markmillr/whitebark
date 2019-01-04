#plotter.py
import matplotlib.pyplot as plt


def plot(x, y, title):

	fig, ax = plt.subplots()
	ax.scatter(x, y)
	ax.plot(x, y)
	ax.set_xlim(0,24)
	ax.set_ylim(-20,400)
	ax.set_xticks(np.arange(0, 25, step=6))
	ax.set_xlabel('Time (lunar hours)')
	ax.set_ylabel('Surface temperature (K)')
	ax.annotate('Max temp: {:.0f} K'.format(np.amax(T_s)), xy=(17,350))
	ax.annotate('Min temp: {:.0f} K'.format(np.amin(T_s)), xy=(17,325))
	ax.set_title(title)

	return plt.show()