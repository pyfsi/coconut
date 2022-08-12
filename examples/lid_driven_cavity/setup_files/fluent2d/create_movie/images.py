import imageio


step = 1
number_of_frames = 700
for variable in ('pressure', 'velocity', 'vector_velocity'):
    file_names = [f'../case_timestep{i * step}-{variable}.jpg' for i in range(1, number_of_frames + 1)]
    with imageio.get_writer(f'movie-{variable}.gif', mode='I') as writer:
        for file_name in file_names:
            image = imageio.imread(file_name)
            writer.append_data(image)
