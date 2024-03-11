from subprocess import Popen, PIPE

p = Popen(['sudo', '-S', 'ls'], stdin=PIPE, stderr=PIPE, stdout=PIPE, text=True)
prompt = p.communicate("password" + '\n')
output = prompt[0]

print(output)