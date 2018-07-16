#!/usr/bin/python3

import tkinter as tk

top = tk.Tk()
top.title('TEST')

top.geometry('600x400')

font = ('Arial', '30')
e1 = tk.Entry(top)
e1.pack()
l1 = tk.Listbox(top)
for i in ('Poa', 'Zea', 'Rosa'):
    l1.insert('end', i)
l1.pack()
r1 = tk.Radiobutton(top, text='option')
r1.pack()
b1 = tk.Button(top, text='Download', font=font)
b2 = tk.Button(top, text='Divide & Filter', font=font)
b3 = tk.Button(top, text='Alignment', font=font)
b4 = tk.Button(top, text='Analysis', font=font)
b5 = tk.Button(top, text='Primer Design', font=font)
b1.pack()
b2.pack()
b3.pack()
b4.pack()
b5.pack()
top.mainloop()
