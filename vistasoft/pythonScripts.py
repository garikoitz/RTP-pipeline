for a in A:
       :     a = a.strip('\n')
       :     extn = pathlib.Path(a).suffix
       :     b = a.replace('toolboxes/vistasoft','soft/vistasoft_garikoitz')
       :     if extn == '.m':
       :         if os.path.isfile(a) and os.path.isfile(b):
       :             print(a)
       :             print(b)
       :             print('---------')
       :             text1 = open(a).readlines()
       :             text2 = open(b).readlines()
       :             for line in difflib.unified_diff(text1, text2):
       :                 print(line,)
       :             print('=========')
       :             print('\n')



# This only copies the .m files, copy the mexa64 files as well
for a in A:
           :     a = a.strip('\n')
       :     extn = pathlib.Path(a).suffix
       :     b = a.replace('toolboxes/vistasoft','soft/vistasoft_garikoitz')
       :     if extn == '.m':
	          :         if os.path.isfile(a) and os.path.isfile(b):
		             :             src = b
       :             dst = a.replace('/data/localhome/glerma/toolboxes/vistasoft','/black/localhome/glerma/soft/RTP-pipeline/vistasoft')
       :             print(src)
       :             print(dst)
       :             if not os.path.isfile(dst): sh.copyfile(src,dst)
       :             print('\n')



