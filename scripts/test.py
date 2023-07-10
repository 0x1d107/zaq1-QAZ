import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from inv_q_filtr import phase_q_correction

signal = np.loadtxt("./scripts/signal.txt")
nperseg = 256
overlap = nperseg // 2

fs = 500

f, t, Zxx = scipy.signal.stft(signal, nperseg=nperseg, noverlap=overlap, fs = 500)

# plt.pcolormesh(t, f, np.abs(Zxx), shading='gouraud')
# plt.show()

plt.plot(signal)
# corrected_signal = phase_q_correction(signal, 500, 200, 20, 256, 128)

# print(len(corrected_signal))
# for i in range(len(corrected_signal)):
#     print(corrected_signal[i])

# plt.plot(signal, label = "aaa")
# # plt.plot(corrected_signal)
# # plt.show()

cor = np.loadtxt("/home/murellos/projects/inv-q-filtr/result.txt")

plt.plot(cor)
# plt.legend()
plt.show()