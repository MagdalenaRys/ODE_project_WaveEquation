import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from io import BytesIO
import Metoda_numeryczna as sm


# Funkcja generująca wykres
def create_plot(values, constant, P, Q):
    m = np.max(np.abs(P[1])) * np.max(np.abs(Q[1])) * 1.2
    x = np.array(P[0])
    y = values * constant

    plt.figure(figsize=(6, 4))
    plt.plot(x, y)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u_{1}(t,x)$')
    plt.title(rf'Drganie struny gitarowej dla $n={n}$')
    plt.grid(True)
    plt.ylim(-m, m)

    image_bytes = BytesIO()
    plt.savefig(image_bytes, format='png')
    plt.close()
    return image_bytes.getvalue()


for n in range(1, 7):
    l = 10.0
    c = 1.0

    d = sm.shootingMethod(lambda w, x, y: sm.funP(w, x, y, n, l), [(0.0, 0.0), (l, 0.0)], (0.0, l), (-3.0, 2.0))
    P = sm.modifiedEulerMethod(lambda w, x, y: sm.funQ(w, x, y, n, l, 1.0), [(0.0, 0.0), (l, c)], (0.0, l))
    Q = sm.modifiedEulerMethod(lambda w, x, y: sm.funQ(w, x, y, n, l, c), sm.initialConditionsToQ(3, sm.f_func, sm.g_func, P),
        (0, 2*l))

    arguments = np.array(P[0])
    values = np.array(P[1])

    a = []
    b = 30*n
    for i in range(b):
        a.append(Q[1][len(Q[1])*i//b])
    constants_list = np.array(a)

    frames = []
    for constant in constants_list:
        frame = create_plot(values, constant, P, Q)
        frames.append(Image.open(BytesIO(frame)))

    # Tworzenie gifa
    gif_bytes = BytesIO()
    frames[0].save(
        gif_bytes,
        format='GIF',
        save_all=True,
        append_images=frames[1:],
        duration=100,
        loop=0
    )

    # Zapisanie gifa
    with open(f'../gifs/struna_gitarowa_n={n}.gif', 'wb') as f:
        f.write(gif_bytes.getvalue())

    print("GIF został zapisany.")
