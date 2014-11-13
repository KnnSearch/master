from Tkinter import *
root = Tk();
c = Canvas(root,width=600,height=600,bg='white');
c.pack();

count = 0
for i in range(50,551,25):
    text = str(count)
    c.create_text(i,35,text = text,justify = CENTER)
    count = count + 1

count = 0
for j in range(50,551,25):
    text = str(count)
    c.create_text(35,j,text = text,justify = CENTER)
    count = count + 1

for i in range(50,551,25):
    ly = c.create_line(i,50,i,550)
    
for j in range(50,551,25):
    lx = c.create_line(50,j,550,j)

data = open('data.txt')
for i in data:
    temp = i[2:-1]
    card_temp = temp.split()
    card = [float(card_temp[0]),float(card_temp[1])]
    point = c.create_oval(card[0]*25+50-5,card[1]*25+50-5,card[0]*25+50+5,card[1]*25+50+5,fill = 'red');

data = open('database_1.txt')
flag = 0;
for i in data:
    if(flag == 0):
        flag = 1
    else:
        temp = i.split()
        temp = temp[1:5]
        rect = [float(temp[0]),float(temp[1]),float(temp[2]),float(temp[3])]
        c.create_rectangle(rect[0]*25+50,rect[1]*25+50,rect[2]*25+50,rect[3]*25+50,width = 5)
'''
data = open('rect.txt')
for i in data:
    temp = i.split()
    rect = [float(temp[0]),float(temp[1]),float(temp[2]),float(temp[3])]
    c.create_rectangle(rect[0]*25+50,rect[1]*25+50,rect[2]*25+50,rect[3]*25+50,width = 5)
'''
input()
