ar=[] # to store the data from the protein.fa file
with open('protein.fa', 'r') as f:
	for line in f:
		ar.append(line[:-1])
# print(ar)
pSeq1=ar[1]
pSeq2=ar[4]
print(pSeq1)
print(pSeq2)
m=len(pSeq1)
n=len(pSeq2)
# print(m,n)
DP=[]
# DP.append((pSeq1))
# 
# Filling the value of 2D array DP with the value 0 
# filling first row with seq 1
# filling first column with seq 2
# so now the size of DP is m x n
for i in range (n+1):
	a=[]
	if(i==0):
		a.append(' ')
		for k in range(m):
			a.append(pSeq1[k])
	else:

		for j in range(1+m):
			if(j==0):
				a.append(pSeq2[i-1])
			else:
				a.append(0)
	DP.append((a))
# print(DP)
# DP initiliazation complete
# 
with open("3rd.txt",'w') as g:

	######################  Dot plot matrix ##############################

	for i in range (n):
		for j in range(m):
			if(pSeq1[j]==pSeq2[i]):
				DP[i+1][j+1]=1

	# printing the dot plot matrix
	for i in range (1+n):
		for j in range(1+m):
			if(DP[i][j]==0):
				print(' ',end=' ')
				print(' ',end=' ',file=g)

			elif(DP[i][j]==1):
				# print(DP[i][j],end=' ')
				print('*',end=' ')
				print('*',end=' ',file=g)

			else:
				print(DP[i][j],end=' ')
				print(DP[i][j],end=' ',file=g)
		print()
	########################Dot plot matrix##########################
	######################  Sum Matrix ##############################

	BTrc=[] # stroing the index(of the max from prev row or column ) from which value is fetched from for the current index 
	# this is 2D array 
	# filling the value at the time of calculating sum matrix
	for i in range (1+n):
		a=[]
		for j in range(1+m):
			a.append([0,0])
		BTrc.append(a)

	for i in range(n,0,-1):
		for j in range(m,0,-1):
			if(i==n or j==m):
				if(pSeq1[j-1]==pSeq2[i-1]):
					DP[i][j]=1
				else:
					DP[i][j]=0
				BTrc[i][j]=[i, j]
			else:
				columnMax=0;rowMax=0;rowMaxIndex=j+1;columnMaxIndex=i+1
				for column in range(j+1,1+m):
					if(DP[i+1][column]>rowMax):
						rowMaxIndex=column
						rowMax=DP[i+1][column]
				for row in range(i+1,1+n):
					if(DP[row][j+1]>columnMax):
						columnMaxIndex=row
						columnMax=DP[row][j+1]
				DP[i][j]+=max(rowMax,columnMax)
				if(rowMax>=columnMax):
					BTrc[i][j][0]=i+1
					BTrc[i][j][1]=rowMaxIndex
				else:
					BTrc[i][j][0]=columnMaxIndex
					BTrc[i][j][1]=j+1

	# printing the sum matrix
	for i in range (1+n):
		for j in range(1+m):
			if(len(str(DP[i][j]))>=2):
				print(DP[i][j],end=' ')
				print(DP[i][j],end=' ',file=g)

			else:
				print(DP[i][j],end='  ')
				print(DP[i][j],end='  ',file=g)


		print()
		print("",file=g)


	######### printing backtracing matrix
	# for i in range (1+n):
	# 	for j in range(1+m):
	# 		print(BTrc[i][j],end='  ')

	# 	print()

	# print(BTrc)
	#################################################################
	maxm=0;maxmIndex=1
	for i in range(1,1+m):
		if(DP[1][i]>maxm):
			maxm=DP[1][i]
			maxmIndex=i
	i=1;j=1;s1='';s2=''
	while(i<n and j<m):
		if(BTrc[i][j][0]==i and BTrc[i][j][1]==j):
			break
		# print(pSeq1[j])
		# print(pSeq2[i])
		s1+=pSeq1[j-1]
		s2+=pSeq2[i-1]
		ind=BTrc[i][j][0] # stroring the value of index using which the current i,jth value is calculated
		jnd=BTrc[i][j][1] # stroring the value of index using which the current i,jth value is calculated
		print("for DP",i,j,"value fetched from",ind,jnd)
		# print("for DP",i,j,"value fetched from",ind,jnd,file=g)
		if(ind!=i+1): #adding the dash '_' if skipping any seq2 char
			for k in range(i+1,ind):
				s2+=pSeq2[k-1]
				s1+="_"
		if(jnd!=j+1): #adding the dash '_' if skipping any seq1 char
			for k in range(j+1,jnd):
				s1+=pSeq1[k-1]
				s2+="_"
		i=ind;j=jnd
	# printing the remaining sequence (either first or second)
	while(i<n):
		s2+=pSeq2[i]
		i+=1
	while(j<m):
		s1+=pSeq1[j]
		j+=1 
	print(s1)
	print(s2)
	print(s1,file=g)
	print(s2,file=g)