import imageio


step = 50
number_of_frames = 350
for variable in ('pressure', 'velocity'):
    file_names = [f'../case_timestep{i * 50}-{variable}.jpg' for i in range(1, number_of_frames + 1)]
    with imageio.get_writer(f'movie-{variable}.gif', mode='I') as writer:
        for file_name in file_names:
            image = imageio.imread(file_name)
            writer.append_data(image)
