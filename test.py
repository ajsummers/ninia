
import re

test_string = '''
   Ag	1.49199531	0.86140389	20.00000000 0 0 0
   Ag	4.47598592	0.86140389	20.00000000 0 0 0
   Ag	2.98399062	3.44561557	20.00000000 0 0 0
   Ag	5.96798123	3.44561557	20.00000000 0 0 0
   Ag	0.00000000	1.72280779	22.43641814 0 0 0
   Ag	2.98399062	1.72280779	22.43641814 0 0 0
   Ag	1.49199531	4.30701946	22.43641814 0 0 0
   Ag	4.47598592	4.30701946	22.43641814 0 0 0
   Ag	0.00000000	0.00000000	24.87283627
   Ag	2.98399062	0.00000000	24.87283627
   Ag	1.49199531	2.58421168	24.87283627
   Ag	4.47598592	2.58421168	24.87283627
   C	1.49199531	0.86140389	26.87283627
'''
lock = [1, 2, 3]
which = (1, 0, 0)
new_string = []
for index, line in enumerate(test_string.splitlines()):
    if index in lock:
        line = re.sub(r' \d \d \d', '', line)
        for value in which:
            line += f' {value}'

    new_string.append(line)

new_string = '\n'.join(new_string)
print(new_string)