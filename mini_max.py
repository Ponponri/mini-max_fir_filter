import argparse
import numpy as np
import math
#from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

class mini_max:
	def __init__(self, args):
		self.filter_length = args.filter_length
		self.fs = args.fs 
		self.pass_band = [args.pass_band_L, args.pass_band_H]
		self.transition_band = [args.transition_band_L,args.transition_band_H]
		self.WFP = args.WFP
		self.WFS = args.WFS
		self.threshold = args.threshold

		self.pass_band_normalize = [x / self.fs for x in self.pass_band]
		self.transition_band_normalize = [x / self.fs for x in self.transition_band]

		self.k = None
		self.mat_A = None
		self.vec_s = None
		self.vec_B = None
		self.Fm = None

	#step 1
	def find_init_Fm(self):
		#Set initial values of Fm
		self.k = (self.filter_length - 1) // 2
		self.mat_A = np.zeros((self.k+2, self.k+2))
		self.vec_s = np.zeros((self.k+2,))
		self.vec_B = np.zeros((self.k+2,))

		if(self.k % 2 == 0):
			FmL = np.linspace(0, self.transition_band_normalize[0], (self.k+2)//2)
			FmH = np.linspace(self.transition_band_normalize[1], 0.5, (self.k+2)//2)
			self.Fm = np.concatenate([FmL, FmH], axis = 0)

		elif(self.k % 2 == 1):
			FmL = np.linspace(0, self.transition_band_normalize[0], (self.k+2)//2)
			FmH = np.linspace(self.transition_band_normalize[1], 0.5, (self.k+2)//2 + 1)
			self.Fm = np.concatenate([FmL, FmH], axis = 0)


	#step 2
	def find_s(self):
		#Set vec Hd
		Hd = np.zeros(self.k+2,)
		for i in range(self.k+2):
			if(self.Fm[i] >= self.pass_band_normalize[0] and self.Fm[i] <= self.pass_band_normalize[1]):
				Hd[i] = 1
			else:
				Hd[i] = 0

		#Set Wf
		Wf = np.zeros(self.k+2,)
		for i in range(self.k+2):
			if(Hd[i] == 1):
				Wf[i] = self.WFP
			else:
				Wf[i] = self.WFS

		#Set mat A
		for i in range(self.k+2):
			for j in range(self.k+2):
				if(j == 0):
					self.mat_A[i][j] = 1
				elif(j == self.k+1):
					self.mat_A[i][j] = (1 / Wf[i]) * math.pow(-1, i)
				else:
					self.mat_A[i][j] = math.cos(2 * j * math.pi * self.Fm[i])

		#Find s = inv(A) * B
		self.vec_s = np.dot(np.linalg.inv(self.mat_A), Hd)


	#step 3~5
	def find_err(self):
		#Find error 
		sample = np.linspace(0,0.5,8000)
		err_list = []
		for i in sample:
			Rf = 0
			for j in range(self.k+1):
				Rf = Rf + self.vec_s[j] * math.cos(2 * j * math.pi * i)
			#print('Rf:',Rf)
			if(i <= self.transition_band_normalize[0] or i >= self.transition_band_normalize[1]):
				if(i >= self.pass_band_normalize[0] and i <= self.pass_band_normalize[1]):
					err_tmp = (Rf - 1) * self.WFP
				else:
					err_tmp = (Rf - 0) * self.WFS

				err_list.append([i, err_tmp])

		#Find local maximum of error
		err_great = []
		for i in range(len(err_list)):
			if(err_list[i][0] == 0):
				if(err_list[i][1] > err_list[i+1][1] and err_list[i][1] > 0):
					err_great.append(err_list[i])
				elif(err_list[i][1] < err_list[i+1][1] and err_list[i][1] < 0):
					err_great.append(err_list[i])

			elif(err_list[i][0] == 0.5):
				if(err_list[i][1] > err_list[i-1][1] and err_list[i][1] > 0):
					err_great.append(err_list[i])
				elif(err_list[i][1] < err_list[i-1][1] and err_list[i][1] < 0):
					err_great.append(err_list[i])

			elif(err_list[i][0] == self.transition_band_normalize[0] or err_list[i][0] == self.transition_band_normalize[1]):
				err_great.append(err_list[i])

			else:
				if(err_list[i][1] > err_list[i-1][1] and err_list[i][1] > err_list[i+1][1]):
					err_great.append(err_list[i])
				elif(err_list[i][1] < err_list[i-1][1] and err_list[i][1] < err_list[i+1][1]):
					err_great.append(err_list[i])

		#Find max error
		err_max = 0
		for i in range(len(err_great)):
			if(err_max < abs(err_great[i][1])):
				err_max = abs(err_great[i][1])

		#Remove points to let length of err_great to k+2 
		while(len(err_great) > self.k+2):
			if((self.pass_band_normalize[0]+self.pass_band_normalize[1]) / 2 >= 0.25):
				del err_great[0]
			else:
				del err_great[-1]

		#Set new Fm on local maximum
		err_great_inx = []
		for i in range(len(err_great)):
			err_great_inx.append(err_great[i][0])
		self.Fm = np.array(err_great_inx)

		return err_max


	#step 6
	def show_response(self,err_list):
		#Get impulse response
		h = np.zeros(self.k * 2 + 1)
		for i in range(len(self.vec_s) - 1):
			if(i == 0):
				h[self.k] = self.vec_s[0]
			else:
				h[self.k+i] = self.vec_s[i]/2
				h[self.k-i] = self.vec_s[i]/2

		#Show results
		super_title = f'filter_length: {self.filter_length}, fs: {self.fs}, pass_band: [{self.pass_band[0]},{self.pass_band[1]}], transition_band: [{self.transition_band[0]},{self.transition_band[1]}],\nweight_function_pass: {self.WFP}, weight_function_stop: {self.WFS},  threshold: {self.threshold}'
		plt.figure(figsize=(10,4))
		plt.suptitle(super_title,fontsize = 8)

		#Plot frequency response
		r = []
		ideal = []
		sample = np.linspace(0,0.5,8000)
		for i in sample:
			Rf = 0
			for j in range(self.k+1):
				Rf = Rf + self.vec_s[j] * math.cos(2 * j * math.pi * i)
			r.append(Rf)

			if(i >= self.pass_band_normalize[0] and i <= self.pass_band_normalize[1]):
				ideal.append(1)
			else:
				ideal.append(0)

		plt.subplot(131)
		plt.title('frequency response')
		plt.xlabel('normal frequency')
		plt.ylabel('r')
		plt.plot(sample, r, label = 'mini-mix')
		plt.plot(sample, ideal, label = 'ideal')
		plt.legend(loc='upper left',fontsize=5)

		#Plot impluse response
		n = np.linspace(0,self.k * 2,self.k * 2 + 1)
		plt.subplot(132)
		plt.title('impluse response')
		plt.xlabel('n')
		plt.ylabel('h')
		plt.bar(n,h)

		#Plot iteration error_max
		x = np.linspace(1,len(err_list),len(err_list))
		plt.subplot(133)
		plt.title('error_max')
		plt.xlabel('iteration')
		plt.ylabel('error')
		plt.xlim(1, len(err_list))
		plt.plot(x,err_list)

		plt.tight_layout()
		plt.savefig('results.png', dpi=300)
		plt.show()








