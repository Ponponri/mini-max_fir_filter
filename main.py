import argparse
from mini_max import mini_max


def main(args):
	mm = mini_max(args)  #Initialize object
	mm.find_init_Fm()  #Set initial values
	mm.find_s()  #Find s = inv(A) * B

	#Get error_max
	err_list = []
	err_current = mm.find_err()
	err_list.append(err_current)

	#Iterative to find results
	while True:
		err_past = err_current
		mm.find_s()
		err_current = mm.find_err()
		err_list.append(err_current)
		if(err_current - err_past <= mm.threshold and err_current - err_past >= 0):
			break
		err_past = err_current

	#Print error_max of each iteration
	for i in range(len(err_list)):
		print(f'iteration {i+1}, err: {err_list[i]}')
	mm.show_response(err_list)
	print('done')

if __name__ == '__main__': 
	#Get parameters
	parser = argparse.ArgumentParser()
	parser.add_argument('--filter_length', default = 21, type = int)
	parser.add_argument('--fs', default = 8000, type = int)
	parser.add_argument('--pass_band_L', default = 1800, type = int)
	parser.add_argument('--pass_band_H', default = 4000, type = int)
	parser.add_argument('--transition_band_L', default = 1600, type = int)
	parser.add_argument('--transition_band_H', default = 2000, type = int)
	parser.add_argument('--WFP', default = 1, type = float)
	parser.add_argument('--WFS', default = 0.8, type = float)
	parser.add_argument('--threshold', default = 0.0001, type = float)

	args = parser.parse_args()
	main(args)