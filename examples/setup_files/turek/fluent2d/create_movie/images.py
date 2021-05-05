import imageio
import os


images_dir = ''
file_names = [f'case_timestep{i + 1}.jpeg' for i in range(len(os.listdir(images_dir)))]
with imageio.get_writer('movie.gif', mode='I') as writer:
    for file_name in file_names:
        image = imageio.imread(os.path.join(images_dir, file_name))
        writer.append_data(image)
