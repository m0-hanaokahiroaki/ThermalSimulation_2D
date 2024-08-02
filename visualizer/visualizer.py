import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches  

#同じディレクトリ内にoutput_xxx.txtファイルを置くことで実行可能

num = 2295#総ファイル数
dt = 0.1#時間ステップ幅
data = np.empty((num, 101, 101), dtype=float)

for i in range(num):
  file_name = f"output_{i}.txt"
  data[i] = np.loadtxt(file_name)


fig, ax = plt.subplots()
image = ax.imshow(data[0], cmap='coolwarm')
plt.colorbar(image)
def animate(i):
  image.set_array(data[i])
  ax.set_title(f'Time {i*dt:.1f}s')
  if i%100 == 0 :
      print(f'{100.0*i/num:.1f}% completed')
  ax.set_xlabel('2x [cm]')
  ax.set_ylabel('2y [cm]')
  return image,

ani = animation.FuncAnimation(fig, animate, frames=num, interval=20, blit=True)
ani.save('animation.gif', writer='pillow')