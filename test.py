from subprocess import Popen, PIPE

p = Popen(['sudo', '-S', 'ls'], stdin=PIPE, stderr=PIPE, stdout=PIPE, text=True)
prompt = p.stdin.write("password" + '\n')
p.stdin.flush()

print(prompt)