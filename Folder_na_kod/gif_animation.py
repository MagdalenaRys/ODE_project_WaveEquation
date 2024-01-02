import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageSequence
from io import BytesIO
import Shooting_Method as sm
import math

P = sm.shootingMethod(lambda w, x, y: sm.funP(w, x, y, math.pi / 1), [(0.0, 0.0), (10.0, 0.0)], (0.0, 10.0), (-3.0, 2.0))

# p = modifiedEulerMethod(lambda w, x, y: funP(w, x, y, math.pi), [(0.0, 0.0), (0.0, 1.0)], (0.0, 10.0))

Q = sm.modifiedEulerMethod(lambda w, x, y: sm.funQ(w, x, y, math.pi / 1, 1.0), sm.initialConditionsToQ(3, sm.f_func, sm.g_func, P),
                        (0.0, 10.0))


# Function to create a plot and return the image as bytes
def create_plot(values, constant):
    x = np.array(P[0])
    y = values * constant

    plt.figure(figsize=(6, 4))
    plt.plot(x, y)
    plt.xlabel('Index')
    plt.ylabel('Values * Constant')
    plt.grid(True)
    plt.tight_layout()

    # Set y-axis limits to maintain consistency
    plt.ylim(-10000, 10000)  # Adjust the limits based on your data range

    # Save the plot as an image in memory
    image_bytes = BytesIO()
    plt.savefig(image_bytes, format='png')
    plt.close()  # Close the plot to avoid displaying multiple plots

    return image_bytes.getvalue()


# Example arrays of arguments and values
arguments = np.array(P[0])
values = np.array(P[1])

# Example list of constants
a = []
n = 30
for i in range(n):
    a.append(Q[1][len(Q[1])*i//n])
constants_list = np.array(a)

# Create frames for each constant
frames = []
for constant in constants_list:
    frame = create_plot(values, constant)
    frames.append(Image.open(BytesIO(frame)))

# Create a GIF from the frames
gif_bytes = BytesIO()
frames[0].save(
    gif_bytes,
    format='GIF',
    save_all=True,
    append_images=frames[1:],
    duration=200,  # Set the duration between frames in milliseconds
    loop=0  # Set loop to 0 for infinite loop, or any other positive integer for a finite loop
)

# Save or display the GIF as needed
# For example, to save the GIF to a file:
with open('multiplication_animation.gif', 'wb') as f:
    f.write(gif_bytes.getvalue())

print("Animated GIF created in memory.")

