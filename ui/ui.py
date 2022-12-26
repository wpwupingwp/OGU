#!/usr/bin/env python
from logging import handlers
from time import time, localtime, strftime
from tkinter import filedialog, messagebox, scrolledtext
import logging
import platform
import queue
import sys
import threading
import tkinter as tk
import tkinter.ttk as ttk
import webbrowser

from BarcodeFinder.global_vars import global_dict, log, FMT, DATEFMT
from BarcodeFinder.gb2fasta import gb2fasta_main
from BarcodeFinder.evaluate import evaluate_main
from BarcodeFinder.primer import primer_main


def set_combo_style(win: tk.Frame):
    style = 'combostyle'
    win.combo_style = ttk.Style()
    if style not in ttk.Style().theme_names():
        win.combo_style.theme_create(style, parent='alt', settings={
            'TCombobox': {'configure': {
                'selectbackground': 'white',
                'fieldbackground': 'white',
                'background': 'white'}}})
    win.combo_style.theme_use(style)


def move_to_center(window: tk.Tk, width: int, height: int) -> None:
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    x = (screen_width - width) // 2
    y = (screen_height - height) // 2
    window.geometry(f'{width}x{height}+{x}+{y}')
    return


def after_close(frame):
    """
    Deiconify root before destroy.
    """

    def func():
        root.deiconify()
        frame.destroy()
    return func


def scroll_text(window):
    """
    ScrolledText that shows logs.
    """
    def poll():
        while True:
            try:
                msg = log_queue.get(block=False)
                level = msg.levelname
                msg = formatter.format(msg) + '\n'
                scroll.insert('end', msg, level)
                scroll.yview('end')
            except queue.Empty:
                break
        # to avoid orphan poll()
        if log.hasHandlers():
            scroll.after(10, poll)
        else:
            return

    # clean old handlers
    for i in log.handlers:
        log.removeHandler(i)
    log_queue = queue.Queue()
    formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)
    # do not add formatter to queuehandler, or msg will be formatted twice
    queue_handler = handlers.QueueHandler(log_queue)
    # give poll() time to quit
    root.after(100, log.addHandler(queue_handler))
    scroll = scrolledtext.ScrolledText(window)
    scroll.tag_config('INFO', foreground='black')
    scroll.tag_config('WARNING', foreground='orange')
    scroll.tag_config('ERROR', foreground='red')
    scroll.tag_config('CRITICAL', foreground='red')
    scroll.tag_config('EXCEPTION', foreground='red')
    scroll.pack(fill='both')
    scroll.after(0, poll)

def thread_wrap(function, arg_str, window):
    """
    Wrap for callback.
    Args:
        function(callable): function to call
        arg_str(str): string for function's argparse
        window(Toplevel): window to hide
    """
    try:
        result = function(arg_str)
    except Exception as e:
        log.exception(str(e))
        log.exception('Abort.')
        messagebox.showinfo(message='Abort.')
        root.deiconify()
        return
    messagebox.showinfo(message=f'Done. See {result[0].out} for details.')
    window.withdraw()
    root.deiconify()
    return


def open_file(entry, single=True, type_='file', entry2=None, title=''):
    """
    Set title, fill entry 1, empty entry 2.
    """

    def func():
        if type_ != 'file':
            a = filedialog.askdirectory(title=title)
        else:
            if single:
                a = filedialog.askopenfilename(title=title)
            else:
                a = filedialog.askopenfilenames(title=title)
        entry.delete(0, 'end')
        entry.insert(0, a)
        if entry2 is not None:
            entry2.delete(0, 'end')

    return func


class Root:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#edf0f3'  # Closest X11 color: 'gray94'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9'  # X11 color: 'gray85'
        _ana1color = '#d9d9d9'  # X11 color: 'gray85'
        _ana2color = '#ececec'  # Closest X11 color: 'gray92'
        _tabfg1 = 'black'
        _tabfg2 = 'black'
        _tabbg1 = 'grey75'
        _tabbg2 = 'grey89'
        _bgmode = 'light'

        photo_location1 = "./button1.png"
        photo_location2 = "./button2.png"
        photo_location3 = "./button3.png"
        photo_location4 = "./button4.png"

        top.geometry("800x450+400+0")
        move_to_center(top, 800, 450)
        top.minsize(120, 15)
        top.resizable(1, 1)
        top.title("BarcodeFinder")
        top.configure(background="#edf0f3")
        top.configure(highlightbackground="#edf0f3")
        top.configure(highlightcolor="black")
        self.top = top

        self.help_b = tk.Button(self.top)
        self.help_b.place(relx=0.913, rely=0.067, height=40, width=40)
        self.help_b.configure(activebackground="#ececec")
        self.help_b.configure(activeforeground="#000000")
        self.help_b.configure(background="#edf0f3")
        self.help_b.configure(borderwidth="1")
        self.help_b.configure(command=run_help)
        self.help_b.configure(foreground="#000000")
        self.help_b.configure(highlightbackground="#edf0f3")
        self.help_b.configure(highlightcolor="black")
        global _img4
        _img4 = tk.PhotoImage(file=photo_location4)
        self.help_b.configure(image=_img4)
        self.help_b.configure(pady="0")
        self.help_b.configure(text='''Button''')

        self.gb2fasta_b = tk.Button(self.top)
        self.gb2fasta_b.place(relx=0.188, rely=0.288, height=100, width=100)
        self.gb2fasta_b.configure(activebackground="#ececec")
        self.gb2fasta_b.configure(activeforeground="#000000")
        self.gb2fasta_b.configure(background="#edf0f3")
        self.gb2fasta_b.configure(borderwidth="0")
        self.gb2fasta_b.configure(command=ui_gb2fasta)
        self.gb2fasta_b.configure(foreground="#000000")
        self.gb2fasta_b.configure(highlightbackground="#edf0f3")
        self.gb2fasta_b.configure(highlightcolor="black")
        global _img0
        _img0 = tk.PhotoImage(file=photo_location1)
        self.gb2fasta_b.configure(image=_img0)
        self.gb2fasta_b.configure(pady="0")
        self.gb2fasta_b.configure(text='''Button''')

        self.evaluate_b = tk.Button(self.top)
        self.evaluate_b.place(relx=0.438, rely=0.288, height=100, width=100)
        self.evaluate_b.configure(activebackground="#ececec")
        self.evaluate_b.configure(activeforeground="#000000")
        self.evaluate_b.configure(background="#edf0f3")
        self.evaluate_b.configure(borderwidth="0")
        self.evaluate_b.configure(command=ui_evaluate)
        self.evaluate_b.configure(foreground="#000000")
        self.evaluate_b.configure(highlightbackground="#edf0f3")
        self.evaluate_b.configure(highlightcolor="black")
        global _img1
        _img1 = tk.PhotoImage(file=photo_location2)
        self.evaluate_b.configure(image=_img1)
        self.evaluate_b.configure(pady="0")
        self.evaluate_b.configure(text='''Button''')

        self.primer_b = tk.Button(self.top)
        self.primer_b.place(relx=0.688, rely=0.288, height=100, width=100)
        self.primer_b.configure(activebackground="#ececec")
        self.primer_b.configure(activeforeground="#000000")
        self.primer_b.configure(background="#edf0f3")
        self.primer_b.configure(borderwidth="0")
        self.primer_b.configure(command=ui_primer)
        self.primer_b.configure(foreground="#000000")
        self.primer_b.configure(highlightbackground="#edf0f3")
        self.primer_b.configure(highlightcolor="black")
        global _img2
        _img2 = tk.PhotoImage(file=photo_location3)
        self.primer_b.configure(image=_img2)
        self.primer_b.configure(pady="0")
        self.primer_b.configure(text='''Button''')

        self.install_third_party = tk.Button(self.top)
        self.install_third_party.place(relx=0.713, rely=0.865, height=30,
                                       width=200)
        self.install_third_party.configure(activebackground="#ececec")
        self.install_third_party.configure(activeforeground="#000000")
        self.install_third_party.configure(background="#edf0f3")
        self.install_third_party.configure(borderwidth="1")
        self.install_third_party.configure(compound='left')
        self.install_third_party.configure(
            font="-family {TkDefaultFont} -size 10")
        self.install_third_party.configure(foreground="#000000")
        self.install_third_party.configure(highlightbackground="#edf0f3")
        self.install_third_party.configure(highlightcolor="black")
        self.install_third_party.configure(pady="0")
        self.install_third_party.configure(
            text='''Install third-party software''')

        self.gb2fasta_label = tk.Label(self.top)
        self.gb2fasta_label.place(relx=0.188, rely=0.532, height=30, width=100)
        self.gb2fasta_label.configure(activebackground="#f9f9f9")
        self.gb2fasta_label.configure(anchor='w')
        self.gb2fasta_label.configure(background="#edf0f3")
        self.gb2fasta_label.configure(compound='left')
        self.gb2fasta_label.configure(font="-family {TKDefaultFont} -size 14")
        self.gb2fasta_label.configure(foreground="#000000")
        self.gb2fasta_label.configure(highlightbackground="#edf0f3")
        self.gb2fasta_label.configure(highlightcolor="black")
        self.gb2fasta_label.configure(text='''GB2Fasta''')

        self.evaluate_label = tk.Label(self.top)
        self.evaluate_label.place(relx=0.45, rely=0.532, height=30, width=100)
        self.evaluate_label.configure(activebackground="#f9f9f9")
        self.evaluate_label.configure(anchor='w')
        self.evaluate_label.configure(background="#edf0f3")
        self.evaluate_label.configure(compound='left')
        self.evaluate_label.configure(font="-family {TKDefaultFont} -size 14")
        self.evaluate_label.configure(foreground="#000000")
        self.evaluate_label.configure(highlightbackground="#edf0f3")
        self.evaluate_label.configure(highlightcolor="black")
        self.evaluate_label.configure(text='''Evaluate''')

        self.primer_label = tk.Label(self.top)
        self.primer_label.place(relx=0.713, rely=0.532, height=30, width=100)
        self.primer_label.configure(activebackground="#f9f9f9")
        self.primer_label.configure(anchor='w')
        self.primer_label.configure(background="#edf0f3")
        self.primer_label.configure(compound='left')
        self.primer_label.configure(font="-family {TKDefaultFont} -size 14")
        self.primer_label.configure(foreground="#000000")
        self.primer_label.configure(highlightbackground="#edf0f3")
        self.primer_label.configure(highlightcolor="black")
        self.primer_label.configure(text='''Primer''')

        self.note_label = tk.Label(self.top)
        self.note_label.place(relx=0.35, rely=0.643, height=35, width=252)
        self.note_label.configure(activebackground="#f9f9f9")
        self.note_label.configure(anchor='w')
        self.note_label.configure(background="#edf0f3")
        self.note_label.configure(compound='left')
        self.note_label.configure(font="-family {TKDefaultFont} -size 14")
        self.note_label.configure(foreground="#000000")
        self.note_label.configure(highlightbackground="#edf0f3")
        self.note_label.configure(highlightcolor="black")
        self.note_label.configure(text='''Click button to run modules''')


class GB2Fasta:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#edf0f3'  # Closest X11 color: 'gray94'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9'  # X11 color: 'gray85'
        _ana1color = '#d9d9d9'  # X11 color: 'gray85'
        _ana2color = '#ececec'  # Closest X11 color: 'gray92'
        _tabfg1 = 'black'
        _tabfg2 = 'black'
        _tabbg1 = 'grey75'
        _tabbg2 = 'grey89'
        _bgmode = 'light'
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.', background=_bgcolor)
        self.style.configure('.', foreground=_fgcolor)
        self.style.configure('.', font="TkDefaultFont")
        self.style.map('.', background=[('selected', _compcolor),
                                        ('active', _ana2color)])
        set_combo_style(self)

        top.geometry("600x800+5+139")
        move_to_center(top, 600, 800)
        top.title("GB2Fasta")
        top.configure(background="#edf0f3")
        top.configure(highlightbackground="#edf0f3")
        top.configure(highlightcolor="black")

        self.top = top
        self.gb = tk.StringVar()
        self.gene = tk.StringVar()
        self.molecular = tk.StringVar()
        self.group = tk.StringVar()
        self.og = tk.StringVar()
        self.refseq = tk.StringVar()
        self.count = tk.StringVar()
        self.min_len = tk.StringVar()
        self.max_len = tk.StringVar()
        self.date_start = tk.StringVar()
        self.date_end = tk.StringVar()
        self.exclude = tk.StringVar()
        self.query = tk.StringVar()
        self.taxon = tk.StringVar()
        self.out = tk.StringVar()
        self.expand = tk.StringVar()
        self.max_name_len = tk.StringVar()
        self.max_gene_len = tk.StringVar()
        self.unique = tk.StringVar()
        self.allow_repeat = tk.IntVar()
        self.allow_invert_repeat = tk.IntVar()
        self.allow_mosaic_repeat = tk.IntVar()
        self.no_divide = tk.IntVar()
        self.rename = tk.IntVar()

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.025, rely=0.013, relheight=0.46
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(font="-family {TkDefaultFont} -size 12")
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Input''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(highlightbackground="#edf0f3")
        self.Labelframe1.configure(highlightcolor="black")
        self.TSeparator1 = ttk.Separator(self.Labelframe1)
        self.TSeparator1.place(relx=0.524, rely=0.19, relheight=0.541
                               , bordermode='ignore')
        self.TSeparator1.configure(orient="vertical")

        self.gbfile_label = tk.Label(self.Labelframe1)
        self.gbfile_label.place(relx=0.050, rely=0.057, height=35, width=100,
                                bordermode='ignore')
        self.gbfile_label.configure(activebackground="#f9f9f9")
        self.gbfile_label.configure(activeforeground="black")
        self.gbfile_label.configure(anchor='w')
        self.gbfile_label.configure(background="#edf0f3")
        self.gbfile_label.configure(compound='left')
        self.gbfile_label.configure(font="-family {TKDefaultFont} -size 12")
        self.gbfile_label.configure(foreground="#000000")
        self.gbfile_label.configure(highlightbackground="#edf0f3")
        self.gbfile_label.configure(highlightcolor="black")
        self.gbfile_label.configure(justify='left')
        self.gbfile_label.configure(text='''Genbank files''')

        self.gb_entry = tk.Entry(self.Labelframe1)
        self.gb_entry.place(relx=0.241, rely=0.057, relheight=0.095,
                            relwidth=0.562, bordermode='ignore')
        self.gb_entry.configure(textvariable=self.gb)
        self.gb_entry.configure(takefocus="")
        self.gb_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.gb_entry_tooltip = ToolTip(self.gb_entry, self.tooltip_font,
                                        '''gb format files''')

        self.gb_file_b = tk.Button(self.Labelframe1)
        self.gb_file_b.place(relx=0.82, rely=0.054, height=35, width=90
                             , bordermode='ignore')
        self.gb_file_b.configure(activebackground="#ececec")
        self.gb_file_b.configure(activeforeground="#000000")
        self.gb_file_b.configure(background="#edf0f3")
        self.gb_file_b.configure(command=open_file(self.gb_entry, single=False))
        self.gb_file_b.configure(compound='left')
        self.gb_file_b.configure(font="-family {TkDefaultFont} -size 12")
        self.gb_file_b.configure(foreground="#000000")
        self.gb_file_b.configure(highlightbackground="#edf0f3")
        self.gb_file_b.configure(highlightcolor="black")
        self.gb_file_b.configure(pady="0")
        self.gb_file_b.configure(text='''Open''')

        self.gene_label = tk.Label(self.Labelframe1)
        self.gene_label.place(relx=0.035, rely=0.217, height=35, width=60,
                              bordermode='ignore')
        self.gene_label.configure(activebackground="#f9f9f9")
        self.gene_label.configure(activeforeground="black")
        self.gene_label.configure(anchor='w')
        self.gene_label.configure(background="#edf0f3")
        self.gene_label.configure(compound='left')
        self.gene_label.configure(font="-family {TkDefaultFont} -size 12")
        self.gene_label.configure(foreground="#000000")
        self.gene_label.configure(highlightbackground="#edf0f3")
        self.gene_label.configure(highlightcolor="black")
        self.gene_label.configure(justify='left')
        self.gene_label.configure(text='''Gene''')

        self.gene_entry = ttk.Entry(self.Labelframe1)
        self.gene_entry.place(relx=0.201, rely=0.217, height=35,
                              relwidth=0.314, bordermode='ignore')
        self.gene_entry.configure(textvariable=self.gene)
        self.gene_entry.configure(takefocus="")
        self.gene_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.gene_entry_tooltip = ToolTip(self.gene_entry, self.tooltip_font,
                                          'gene name')

        self.taxon_label = tk.Label(self.Labelframe1)
        self.taxon_label.place(relx=0.035, rely=0.326, height=35, width=70
                               , bordermode='ignore')
        self.taxon_label.configure(activebackground="#f9f9f9")
        self.taxon_label.configure(activeforeground="black")
        self.taxon_label.configure(anchor='w')
        self.taxon_label.configure(background="#edf0f3")
        self.taxon_label.configure(compound='left')
        self.taxon_label.configure(font="-family {TkDefaultFont} -size 12")
        self.taxon_label.configure(foreground="#000000")
        self.taxon_label.configure(highlightbackground="#edf0f3")
        self.taxon_label.configure(highlightcolor="black")
        self.taxon_label.configure(justify='left')
        self.taxon_label.configure(text='''Taxonomy''')

        self.molecular_label = tk.Label(self.Labelframe1)
        self.molecular_label.place(relx=0.558, rely=0.217, height=35, width=80
                                   , bordermode='ignore')
        self.molecular_label.configure(activebackground="#f9f9f9")
        self.molecular_label.configure(activeforeground="black")
        self.molecular_label.configure(anchor='w')
        self.molecular_label.configure(background="#edf0f3")
        self.molecular_label.configure(compound='left')
        self.molecular_label.configure(font="-family {TkDefaultFont} -size 12")
        self.molecular_label.configure(foreground="#000000")
        self.molecular_label.configure(highlightbackground="#edf0f3")
        self.molecular_label.configure(highlightcolor="black")
        self.molecular_label.configure(justify='left')
        self.molecular_label.configure(text='''Molecular''')

        self.TCombobox_molecular = ttk.Combobox(self.Labelframe1)
        self.TCombobox_molecular.place(relx=0.716, rely=0.217, height=35,
                                       width=150, bordermode='ignore')
        self.value_list = ['all', 'DNA', 'RNA', ]
        self.TCombobox_molecular.configure(values=self.value_list)
        self.TCombobox_molecular.configure(state='readonly')
        self.TCombobox_molecular.configure(textvariable=self.molecular)
        self.TCombobox_molecular.current(0)
        self.TCombobox_molecular.configure(takefocus="")

        self.group_label = tk.Label(self.Labelframe1)
        self.group_label.place(relx=0.558, rely=0.326, height=35, width=60
                               , bordermode='ignore')
        self.group_label.configure(activebackground="#f9f9f9")
        self.group_label.configure(activeforeground="black")
        self.group_label.configure(anchor='w')
        self.group_label.configure(background="#edf0f3")
        self.group_label.configure(compound='left')
        self.group_label.configure(font="-family {TkDefaultFont} -size 12")
        self.group_label.configure(foreground="#000000")
        self.group_label.configure(highlightbackground="#edf0f3")
        self.group_label.configure(highlightcolor="black")
        self.group_label.configure(justify='left')
        self.group_label.configure(text='''Group''')

        self.TCombobox_group = ttk.Combobox(self.Labelframe1)
        self.TCombobox_group.place(relx=0.716, rely=0.326, height=35, width=150,
                                   bordermode='ignore')
        self.value_list = ['all', 'animals', 'plants', 'fungi', 'protists',
                           'bacteria', 'archaea', 'viruses', ]
        self.TCombobox_group.configure(values=self.value_list)
        self.TCombobox_group.configure(state='readonly')
        self.TCombobox_group.configure(textvariable=self.group)
        self.TCombobox_group.configure(takefocus="")
        self.TCombobox_group.current(0)

        self.organelle_label = tk.Label(self.Labelframe1)
        self.organelle_label.place(relx=0.035, rely=0.435, height=35, width=80
                                   , bordermode='ignore')
        self.organelle_label.configure(activebackground="#f9f9f9")
        self.organelle_label.configure(activeforeground="black")
        self.organelle_label.configure(anchor='w')
        self.organelle_label.configure(background="#edf0f3")
        self.organelle_label.configure(compound='left')
        self.organelle_label.configure(font="-family {TkDefaultFont} -size 12")
        self.organelle_label.configure(foreground="#000000")
        self.organelle_label.configure(highlightbackground="#edf0f3")
        self.organelle_label.configure(highlightcolor="black")
        self.organelle_label.configure(justify='left')
        self.organelle_label.configure(text='''Organelle''')

        self.TCombobox_og = ttk.Combobox(self.Labelframe1)
        self.TCombobox_og.place(relx=0.201, rely=0.435, height=35, width=170,
                                bordermode='ignore')
        self.value_list = ['ignore', 'both', 'no', 'mitochondrion',
                           'chloroplast', ]
        self.TCombobox_og.configure(values=self.value_list)
        self.TCombobox_og.configure(state='readonly')
        self.TCombobox_og.configure(textvariable=self.og)
        self.TCombobox_og.configure(takefocus="")
        self.TCombobox_og.current(0)

        self.refseq_label = tk.Label(self.Labelframe1)
        self.refseq_label.place(relx=0.035, rely=0.541, height=35, width=70
                                , bordermode='ignore')
        self.refseq_label.configure(activebackground="#f9f9f9")
        self.refseq_label.configure(activeforeground="black")
        self.refseq_label.configure(anchor='w')
        self.refseq_label.configure(background="#edf0f3")
        self.refseq_label.configure(compound='left')
        self.refseq_label.configure(font="-family {TkDefaultFont} -size 12")
        self.refseq_label.configure(foreground="#000000")
        self.refseq_label.configure(highlightbackground="#edf0f3")
        self.refseq_label.configure(highlightcolor="black")
        self.refseq_label.configure(justify='left')
        self.refseq_label.configure(text='''RefSeq''')

        self.TCombobox_refseq = ttk.Combobox(self.Labelframe1)
        self.TCombobox_refseq.place(relx=0.201, rely=0.541, height=35,
                                    width=150, bordermode='ignore')
        self.value_list = ['both', 'yes', 'no', ]
        self.TCombobox_refseq.configure(values=self.value_list)
        self.TCombobox_refseq.configure(state='readonly')
        self.TCombobox_refseq.configure(textvariable=self.refseq)
        self.TCombobox_refseq.configure(takefocus="")
        self.TCombobox_refseq.current(0)
        self.tooltip_font = "TkDefaultFont"
        self.TCombobox_refseq_tooltip = ToolTip(self.TCombobox_refseq,
                                                self.tooltip_font,
                                                'Use RefSeq or not')

        self.count_label = tk.Label(self.Labelframe1)
        self.count_label.place(relx=0.035, rely=0.649, height=35, width=60
                               , bordermode='ignore')
        self.count_label.configure(activebackground="#f9f9f9")
        self.count_label.configure(activeforeground="black")
        self.count_label.configure(anchor='w')
        self.count_label.configure(background="#edf0f3")
        self.count_label.configure(compound='left')
        self.count_label.configure(font="-family {TkDefaultFont} -size 12")
        self.count_label.configure(foreground="#000000")
        self.count_label.configure(highlightbackground="#edf0f3")
        self.count_label.configure(highlightcolor="black")
        self.count_label.configure(justify='left')
        self.count_label.configure(text='''Count''')

        self.count_entry = ttk.Entry(self.Labelframe1)
        self.count_entry.place(relx=0.201, rely=0.649, relheight=0.095
                               , relwidth=0.314, bordermode='ignore')
        self.count_entry.configure(textvariable=self.count)
        self.count_entry.configure(takefocus="")
        self.count_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.count_entry_tooltip = \
            ToolTip(self.count_entry, self.tooltip_font,
                    '''numbers of sequences to download, 0 for no limit''')

        self.len_label = tk.Label(self.Labelframe1)
        self.len_label.place(relx=0.558, rely=0.435, height=35, width=60
                             , bordermode='ignore')
        self.len_label.configure(activebackground="#f9f9f9")
        self.len_label.configure(activeforeground="black")
        self.len_label.configure(anchor='w')
        self.len_label.configure(background="#edf0f3")
        self.len_label.configure(compound='left')
        self.len_label.configure(font="-family {TkDefaultFont} -size 12")
        self.len_label.configure(foreground="#000000")
        self.len_label.configure(highlightbackground="#edf0f3")
        self.len_label.configure(highlightcolor="black")
        self.len_label.configure(justify='left')
        self.len_label.configure(text='''Length''')

        self.min_len_entry = ttk.Entry(self.Labelframe1)
        self.min_len_entry.place(relx=0.686, rely=0.435, relheight=0.095
                                 , relwidth=0.087, bordermode='ignore')
        self.min_len_entry.configure(textvariable=self.min_len)
        self.min_len_entry.configure(takefocus="")
        self.min_len_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.min_len_entry_tooltip = \
            ToolTip(self.min_len_entry, self.tooltip_font,
                    '''sequence length limit''')

        self.to_label = tk.Label(self.Labelframe1)
        self.to_label.place(relx=0.785, rely=0.435, height=35, width=36
                            , bordermode='ignore')
        self.to_label.configure(activebackground="#f9f9f9")
        self.to_label.configure(activeforeground="SystemButtonText")
        self.to_label.configure(anchor='w')
        self.to_label.configure(background="#edf0f3")
        self.to_label.configure(compound='left')
        self.to_label.configure(font="-family {TkDefaultFont} -size 13")
        self.to_label.configure(foreground="#000000")
        self.to_label.configure(highlightbackground="#edf0f3")
        self.to_label.configure(highlightcolor="black")
        self.to_label.configure(text='''to''')

        self.max_len_entry = ttk.Entry(self.Labelframe1)
        self.max_len_entry.place(relx=0.839, rely=0.435, relheight=0.095
                                 , relwidth=0.14, bordermode='ignore')
        self.max_len_entry.configure(textvariable=self.max_len)
        self.max_len_entry.configure(takefocus="")
        self.max_len_entry.configure(cursor="fleur")
        self.min_len_entry.insert(0, '0')
        self.max_len_entry.insert(0, '300000')

        self.date_label = tk.Label(self.Labelframe1)
        self.date_label.place(relx=0.558, rely=0.541, height=35, width=60
                              , bordermode='ignore')
        self.date_label.configure(activebackground="#f9f9f9")
        self.date_label.configure(activeforeground="black")
        self.date_label.configure(anchor='w')
        self.date_label.configure(background="#edf0f3")
        self.date_label.configure(compound='left')
        self.date_label.configure(font="-family {TkDefaultFont} -size 12")
        self.date_label.configure(foreground="#000000")
        self.date_label.configure(highlightbackground="#edf0f3")
        self.date_label.configure(highlightcolor="black")
        self.date_label.configure(justify='left')
        self.date_label.configure(text='''Date''')

        self.date_start_entry = ttk.Entry(self.Labelframe1)
        self.date_start_entry.place(relx=0.686, rely=0.541, relheight=0.095
                                    , relwidth=0.087, bordermode='ignore')
        self.date_start_entry.configure(textvariable=self.date_start)
        self.date_start_entry.configure(takefocus="")
        self.date_start_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.date_start_entry_tooltip = \
            ToolTip(self.date_start_entry, self.tooltip_font, '''1970/1/1''')

        self.to2_label = tk.Label(self.Labelframe1)
        self.to2_label.place(relx=0.789, rely=0.541, height=35, width=36
                             , bordermode='ignore')
        self.to2_label.configure(activebackground="#f9f9f9")
        self.to2_label.configure(activeforeground="SystemButtonText")
        self.to2_label.configure(anchor='w')
        self.to2_label.configure(background="#edf0f3")
        self.to2_label.configure(compound='left')
        self.to2_label.configure(font="-family {TkDefaultFont} -size 13")
        self.to2_label.configure(foreground="#000000")
        self.to2_label.configure(highlightbackground="#edf0f3")
        self.to2_label.configure(highlightcolor="black")
        self.to2_label.configure(text='''to''')

        self.date_end_entry = ttk.Entry(self.Labelframe1)
        self.date_end_entry.place(relx=0.839, rely=0.541, relheight=0.095
                                  , relwidth=0.136, bordermode='ignore')
        self.date_end_entry.configure(textvariable=self.date_end)
        self.date_end_entry.configure(takefocus="")
        self.date_end_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.date_end_entry_tooltip = \
            ToolTip(self.date_end_entry, self.tooltip_font, '''2022/12/31''')

        self.exclude_label = tk.Label(self.Labelframe1)
        self.exclude_label.place(relx=0.558, rely=0.649, height=35, width=60
                                 , bordermode='ignore')
        self.exclude_label.configure(activebackground="#f9f9f9")
        self.exclude_label.configure(activeforeground="black")
        self.exclude_label.configure(anchor='w')
        self.exclude_label.configure(background="#edf0f3")
        self.exclude_label.configure(compound='left')
        self.exclude_label.configure(font="-family {TkDefaultFont} -size 12")
        self.exclude_label.configure(foreground="#000000")
        self.exclude_label.configure(highlightbackground="#edf0f3")
        self.exclude_label.configure(highlightcolor="black")
        self.exclude_label.configure(justify='left')
        self.exclude_label.configure(text='''Exclude''')

        self.exclude_entry = ttk.Entry(self.Labelframe1)
        self.exclude_entry.place(relx=0.686, rely=0.649, relheight=0.095
                                 , relwidth=0.279, bordermode='ignore')
        self.exclude_entry.configure(textvariable=self.exclude)
        self.exclude_entry.configure(takefocus="")
        self.exclude_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.exclude_entry_tooltip = \
            ToolTip(self.exclude_entry, self.tooltip_font,
                    '''exclude expression''')

        self.query_label = tk.Label(self.Labelframe1)
        self.query_label.place(relx=0.035, rely=0.813, height=35, width=60
                               , bordermode='ignore')
        self.query_label.configure(activebackground="#f9f9f9")
        self.query_label.configure(activeforeground="black")
        self.query_label.configure(anchor='w')
        self.query_label.configure(background="#edf0f3")
        self.query_label.configure(compound='left')
        self.query_label.configure(font="-family {TkDefaultFont} -size 12")
        self.query_label.configure(foreground="#000000")
        self.query_label.configure(highlightbackground="#edf0f3")
        self.query_label.configure(highlightcolor="black")
        self.query_label.configure(justify='left')
        self.query_label.configure(text='''Query''')

        self.query_entry = ttk.Entry(self.Labelframe1)
        self.query_entry.place(relx=0.201, rely=0.813, relheight=0.095
                               , relwidth=0.785, bordermode='ignore')
        self.query_entry.configure(textvariable=self.query)
        self.query_entry.configure(takefocus="")
        self.query_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.query_entry_tooltip = ToolTip(self.query_entry, self.tooltip_font,
                                           '''Entrez query string''')

        self.taxon_entry = ttk.Entry(self.Labelframe1)
        self.taxon_entry.place(relx=0.201, rely=0.326, relheight=0.095
                               , relwidth=0.314, bordermode='ignore')
        self.taxon_entry.configure(textvariable=self.taxon)
        self.taxon_entry.configure(takefocus="")
        self.taxon_entry.configure(cursor="fleur")

        self.out_label = tk.Label(self.top)
        self.out_label.place(relx=0.05, rely=0.563, height=36
                             , width=60)
        self.out_label.configure(activebackground="#f9f9f9")
        self.out_label.configure(activeforeground="black")
        self.out_label.configure(anchor='w')
        self.out_label.configure(background="#edf0f3")
        self.out_label.configure(compound='left')
        self.out_label.configure(font="-family {TkDefaultFont} -size 12")
        self.out_label.configure(foreground="#000000")
        self.out_label.configure(highlightbackground="#edf0f3")
        self.out_label.configure(highlightcolor="black")
        self.out_label.configure(justify='left')
        self.out_label.configure(text='''Output''')

        self.out_entry = ttk.Entry(self.top)
        self.out_entry.place(relx=0.167, rely=0.563, relheight=0.045
                             , relwidth=0.633)
        self.out_entry.configure(textvariable=self.out)
        self.out_entry.configure(takefocus="")
        self.out_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.out_entry_tooltip = ToolTip(self.out_entry, self.tooltip_font,
                                         'output folder')
        self.out_b = tk.Button(self.top)
        self.out_b.place(relx=0.833, rely=0.563, height=35, width=80)
        self.out_b.configure(activebackground="#edf0f3")
        self.out_b.configure(activeforeground="#000000")
        self.out_b.configure(background="#edf0f3")
        self.out_b.configure(command=open_file(self.out_entry, type_='folder'))
        self.out_b.configure(compound='left')
        self.out_b.configure(font="-family {TkDefaultFont} -size 12")
        self.out_b.configure(foreground="#000000")
        self.out_b.configure(highlightbackground="#edf0f3")
        self.out_b.configure(highlightcolor="black")
        self.out_b.configure(pady="0")
        self.out_b.configure(text='''Open''')

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.025, rely=0.65, relheight=0.176
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(font="-family {TkDefaultFont} -size 14")
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Advance''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(highlightbackground="#edf0f3")
        self.Labelframe1.configure(highlightcolor="black")

        self.allow_repeat_button = tk.Checkbutton(self.Labelframe1)
        self.allow_repeat_button.place(relx=0.576, rely=0.142, relheight=0.248
                                       , relwidth=0.276, bordermode='ignore')
        self.allow_repeat_button.configure(activebackground="#ececec")
        self.allow_repeat_button.configure(activeforeground="#000000")
        self.allow_repeat_button.configure(anchor='w')
        self.allow_repeat_button.configure(background="#edf0f3")
        self.allow_repeat_button.configure(compound='left')
        self.allow_repeat_button.configure(
            font="-family {TkDefaultFont} -size 12")
        self.allow_repeat_button.configure(foreground="#000000")
        self.allow_repeat_button.configure(highlightbackground="#edf0f3")
        self.allow_repeat_button.configure(highlightcolor="black")
        self.allow_repeat_button.configure(justify='left')
        self.allow_repeat_button.configure(selectcolor="#d9d9d9")
        self.allow_repeat_button.configure(text='''Allow repeat''')
        self.allow_repeat_button.configure(variable=self.allow_repeat)

        self.allow_invert_repeat_button = tk.Checkbutton(self.Labelframe1)
        self.allow_invert_repeat_button.place(relx=0.576, rely=0.709,
                                              relheight=0.248
                                              , relwidth=0.361,
                                              bordermode='ignore')
        self.allow_invert_repeat_button.configure(activebackground="#ececec")
        self.allow_invert_repeat_button.configure(activeforeground="#000000")
        self.allow_invert_repeat_button.configure(anchor='w')
        self.allow_invert_repeat_button.configure(background="#edf0f3")
        self.allow_invert_repeat_button.configure(compound='left')
        self.allow_invert_repeat_button.configure(
            font="-family {TkDefaultFont} -size 12")
        self.allow_invert_repeat_button.configure(foreground="#000000")
        self.allow_invert_repeat_button.configure(highlightbackground="#edf0f3")
        self.allow_invert_repeat_button.configure(highlightcolor="black")
        self.allow_invert_repeat_button.configure(justify='left')
        self.allow_invert_repeat_button.configure(selectcolor="#d9d9d9")
        self.allow_invert_repeat_button.configure(
            text='''Allow invert repeat''')
        self.allow_invert_repeat_button.configure(
            variable=self.allow_invert_repeat)

        self.allow_mosaic_button = tk.Checkbutton(self.Labelframe1)
        self.allow_mosaic_button.place(relx=0.576, rely=0.426, relheight=0.248
                                       , relwidth=0.361, bordermode='ignore')
        self.allow_mosaic_button.configure(activebackground="#ececec")
        self.allow_mosaic_button.configure(activeforeground="#000000")
        self.allow_mosaic_button.configure(anchor='w')
        self.allow_mosaic_button.configure(background="#edf0f3")
        self.allow_mosaic_button.configure(compound='left')
        self.allow_mosaic_button.configure(
            font="-family {TkDefaultFont} -size 12")
        self.allow_mosaic_button.configure(foreground="#000000")
        self.allow_mosaic_button.configure(highlightbackground="#edf0f3")
        self.allow_mosaic_button.configure(highlightcolor="black")
        self.allow_mosaic_button.configure(justify='left')
        self.allow_mosaic_button.configure(selectcolor="#d9d9d9")
        self.allow_mosaic_button.configure(text='''Allow mosaic repeat''')
        self.allow_mosaic_button.configure(variable=self.allow_mosaic_repeat)

        self.expand_label = tk.Label(self.Labelframe1)
        self.expand_label.place(relx=0.175, rely=0.142, height=35
                                , width=60, bordermode='ignore')
        self.expand_label.configure(activebackground="#f9f9f9")
        self.expand_label.configure(activeforeground="black")
        self.expand_label.configure(anchor='e')
        self.expand_label.configure(background="#edf0f3")
        self.expand_label.configure(compound='left')
        self.expand_label.configure(font="-family {TkDefaultFont} -size 12")
        self.expand_label.configure(foreground="#000000")
        self.expand_label.configure(highlightbackground="#edf0f3")
        self.expand_label.configure(highlightcolor="black")
        self.expand_label.configure(justify='left')
        self.expand_label.configure(text='''Expand''')

        self.expand_entry = ttk.Entry(self.Labelframe1)
        self.expand_entry.place(relx=0.297, rely=0.142, relheight=0.248
                                , relwidth=0.187, bordermode='ignore')
        self.expand_entry.configure(textvariable=self.expand)
        self.expand_entry.configure(takefocus="")
        self.expand_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.expand_entry_tooltip = \
            ToolTip(self.expand_entry, self.tooltip_font,
                    '''expand for primer design''')

        self.max_name_len_label = tk.Label(self.Labelframe1)
        self.max_name_len_label.place(relx=0.05, rely=0.426, height=35
                                      , width=130, bordermode='ignore')
        self.max_name_len_label.configure(activebackground="#f9f9f9")
        self.max_name_len_label.configure(activeforeground="black")
        self.max_name_len_label.configure(anchor='e')
        self.max_name_len_label.configure(background="#edf0f3")
        self.max_name_len_label.configure(compound='left')
        self.max_name_len_label.configure(
            font="-family {TkDefaultFont} -size 12")
        self.max_name_len_label.configure(foreground="#000000")
        self.max_name_len_label.configure(highlightbackground="#edf0f3")
        self.max_name_len_label.configure(highlightcolor="black")
        self.max_name_len_label.configure(justify='left')
        self.max_name_len_label.configure(text='''Max name length''')

        self.max_frag_len_label = tk.Label(self.Labelframe1)
        self.max_frag_len_label.place(relx=0.052, rely=0.709, height=35
                                      , width=130, bordermode='ignore')
        self.max_frag_len_label.configure(activebackground="#f9f9f9")
        self.max_frag_len_label.configure(activeforeground="black")
        self.max_frag_len_label.configure(anchor='e')
        self.max_frag_len_label.configure(background="#edf0f3")
        self.max_frag_len_label.configure(compound='left')
        self.max_frag_len_label.configure(
            font="-family {TkDefaultFont} -size 12")
        self.max_frag_len_label.configure(foreground="#000000")
        self.max_frag_len_label.configure(highlightbackground="#edf0f3")
        self.max_frag_len_label.configure(highlightcolor="black")
        self.max_frag_len_label.configure(justify='left')
        self.max_frag_len_label.configure(text='''Max fragment length''')

        self.max_name_len_entry = ttk.Entry(self.Labelframe1)
        self.max_name_len_entry.place(relx=0.297, rely=0.426, relheight=0.248
                                      , relwidth=0.187, bordermode='ignore')
        self.max_name_len_entry.configure(textvariable=self.max_name_len)
        self.max_name_len_entry.configure(takefocus="")
        self.max_name_len_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.max_name_len_entry_tooltip = \
            ToolTip(self.max_name_len_entry, self.tooltip_font,
                    '''max feature name length''')

        self.max_gene_len_entry = ttk.Entry(self.Labelframe1)
        self.max_gene_len_entry.place(relx=0.297, rely=0.709, relheight=0.248
                                      , relwidth=0.187, bordermode='ignore')
        self.max_gene_len_entry.configure(textvariable=self.max_gene_len)
        self.max_gene_len_entry.configure(takefocus="")
        self.max_gene_len_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.max_gene_len_entry_tooltip = \
            ToolTip(self.max_gene_len_entry, self.tooltip_font,
                    '''max fragment sequence length''')
        self.max_name_len_entry.insert(0, '100')
        self.max_gene_len_entry.insert(0, '20000')

        self.TSeparator2 = ttk.Separator(self.Labelframe1)
        self.TSeparator2.place(relx=0.524, rely=0.135, relheight=0.78
                               , bordermode='ignore')
        self.TSeparator2.configure(orient="vertical")

        self.no_divide_b = tk.Checkbutton(self.top)
        self.no_divide_b.place(relx=0.733, rely=0.488, relheight=0.044
                               , relwidth=0.215)
        self.no_divide_b.configure(activebackground="#ececec")
        self.no_divide_b.configure(activeforeground="#000000")
        self.no_divide_b.configure(anchor='w')
        self.no_divide_b.configure(background="#edf0f3")
        self.no_divide_b.configure(compound='left')
        self.no_divide_b.configure(font="-family {TkDefaultFont} -size 12")
        self.no_divide_b.configure(foreground="#000000")
        self.no_divide_b.configure(highlightbackground="#edf0f3")
        self.no_divide_b.configure(highlightcolor="black")
        self.no_divide_b.configure(justify='left')
        self.no_divide_b.configure(selectcolor="#d9d9d9")
        self.no_divide_b.configure(text='''No divide''')
        self.no_divide_b.configure(variable=self.no_divide)

        self.rename_b = tk.Checkbutton(self.top)
        self.rename_b.place(relx=0.467, rely=0.488, relheight=0.044
                            , relwidth=0.243)
        self.rename_b.configure(activebackground="#ececec")
        self.rename_b.configure(activeforeground="#000000")
        self.rename_b.configure(anchor='w')
        self.rename_b.configure(background="#edf0f3")
        self.rename_b.configure(compound='left')
        self.rename_b.configure(font="-family {TkDefaultFont} -size 12")
        self.rename_b.configure(foreground="#000000")
        self.rename_b.configure(highlightbackground="#edf0f3")
        self.rename_b.configure(highlightcolor="black")
        self.rename_b.configure(justify='left')
        self.rename_b.configure(selectcolor="#d9d9d9")
        self.rename_b.configure(text='''Rename gene''')
        self.rename_b.configure(variable=self.rename)

        self.unique_label = tk.Label(self.top)
        self.unique_label.place(relx=0.05, rely=0.488, height=35
                                , width=69)
        self.unique_label.configure(activebackground="#f9f9f9")
        self.unique_label.configure(activeforeground="black")
        self.unique_label.configure(anchor='w')
        self.unique_label.configure(background="#edf0f3")
        self.unique_label.configure(compound='left')
        self.unique_label.configure(font="-family {TkDefaultFont} -size 12")
        self.unique_label.configure(foreground="#000000")
        self.unique_label.configure(highlightbackground="#edf0f3")
        self.unique_label.configure(highlightcolor="black")
        self.unique_label.configure(justify='left')
        self.unique_label.configure(text='''Unique''')

        self.TCombobox_unique = ttk.Combobox(self.top)
        self.TCombobox_unique.place(relx=0.167, rely=0.488, relheight=0.044
                                    , relwidth=0.245)
        self.value_list = ['longest', 'first', 'no', ]
        self.TCombobox_unique.configure(values=self.value_list)
        self.TCombobox_unique.configure(state='readonly')
        self.TCombobox_unique.configure(textvariable=self.unique)
        self.TCombobox_unique.configure(takefocus="")
        self.TCombobox_unique.current(0)
        self.TCombobox_unique_tooltip = ToolTip(
            self.TCombobox_unique, self.tooltip_font,
            'methods to remove redundant records')

        self.run_b = tk.Button(self.top)
        self.run_b.place(relx=0.333, rely=0.863, height=40, width=189)
        self.run_b.configure(activebackground="#ececec")
        self.run_b.configure(activeforeground="#000000")
        self.run_b.configure(background="#edf0f3")
        self.run_b.configure(command=run_gb2fasta(self, self.top))
        self.run_b.configure(compound='left')
        self.run_b.configure(font="-family {TkDefaultFont} -size 14")
        self.run_b.configure(foreground="#000000")
        self.run_b.configure(highlightbackground="#edf0f3")
        self.run_b.configure(highlightcolor="black")
        self.run_b.configure(pady="0")
        self.run_b.configure(text='''Run''')


class Evaluate:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#edf0f3'  # Closest X11 color: 'gray94'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9'  # X11 color: 'gray85'
        _ana1color = '#d9d9d9'  # X11 color: 'gray85'
        _ana2color = '#ececec'  # Closest X11 color: 'gray92'
        _tabfg1 = 'black'
        _tabfg2 = 'black'
        _tabbg1 = 'grey75'
        _tabbg2 = 'grey89'
        _bgmode = 'light'
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.', background=_bgcolor)
        self.style.configure('.', foreground=_fgcolor)
        self.style.configure('.', font="TkDefaultFont")
        self.style.map('.', background=[('selected', _compcolor),
                                        ('active', _ana2color)])

        top.geometry("600x450+109+248")
        move_to_center(top, 600, 450)
        top.title("Evaluate")
        top.configure(background="#edf0f3")
        top.configure(highlightbackground="#d9d9d9")
        top.configure(highlightcolor="black")

        self.top = top
        self.fasta = tk.StringVar()
        self.fasta_folder = tk.StringVar()
        self.aln = tk.StringVar()
        self.out = tk.StringVar()
        self.size = tk.StringVar()
        self.step = tk.StringVar()
        self.quick = tk.IntVar()
        self.ig = tk.IntVar()
        self.iab = tk.IntVar()

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.0, relheight=0.333
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Input''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(highlightbackground="#d9d9d9")
        self.Labelframe1.configure(highlightcolor="black")

        self.unalign_label = tk.Label(self.Labelframe1)
        self.unalign_label.place(relx=0.03, rely=0.133, height=35, width=160
                                 , bordermode='ignore')
        self.unalign_label.configure(activebackground="#f9f9f9")
        self.unalign_label.configure(activeforeground="black")
        self.unalign_label.configure(anchor='w')
        self.unalign_label.configure(background="#edf0f3")
        self.unalign_label.configure(compound='left')
        self.unalign_label.configure(font="-family {TkDefaultFont} -size 10")
        self.unalign_label.configure(foreground="#000000")
        self.unalign_label.configure(highlightbackground="#edf0f3")
        self.unalign_label.configure(highlightcolor="black")
        self.unalign_label.configure(justify='left')
        self.unalign_label.configure(text='''Unaligned FASTA files''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_tooltip = \
            ToolTip(self.unalign_label, self.tooltip_font, '''unaligned''')

        self.fasta_entry = ttk.Entry(self.Labelframe1)
        self.fasta_entry.place(relx=0.314, rely=0.133, relheight=0.233
                               , relwidth=0.489, bordermode='ignore')
        self.fasta_entry.configure(textvariable=self.fasta)
        self.fasta_entry.configure(takefocus="")
        self.fasta_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.fasta_entry_tooltip = \
            ToolTip(self.fasta_entry, self.tooltip_font,
                    '''unaligned fasta files''')

        self.open_btn = tk.Button(self.Labelframe1)
        self.open_btn.place(relx=0.82, rely=0.133, height=35, width=90
                            , bordermode='ignore')
        self.open_btn.configure(activebackground="#ececec")
        self.open_btn.configure(activeforeground="#000000")
        self.open_btn.configure(background="#edf0f3")
        self.open_btn.configure(command=open_file(self.fasta_entry,
                                                  single=False))
        self.open_btn.configure(compound='left')
        self.open_btn.configure(font="-family {TkDefaultFont} -size 10")
        self.open_btn.configure(foreground="#000000")
        self.open_btn.configure(highlightbackground="#edf0f3")
        self.open_btn.configure(highlightcolor="black")
        self.open_btn.configure(pady="0")
        self.open_btn.configure(text='''Open''')

        self.unalign_label2 = tk.Label(self.Labelframe1)
        self.unalign_label2.place(relx=0.03, rely=0.4, height=35, width=160
                                  , bordermode='ignore')
        self.unalign_label2.configure(activebackground="#f9f9f9")
        self.unalign_label2.configure(activeforeground="black")
        self.unalign_label2.configure(anchor='w')
        self.unalign_label2.configure(background="#edf0f3")
        self.unalign_label2.configure(compound='left')
        self.unalign_label2.configure(font="-family {TkDefaultFont} -size 10")
        self.unalign_label2.configure(foreground="#000000")
        self.unalign_label2.configure(highlightbackground="#edf0f3")
        self.unalign_label2.configure(highlightcolor="black")
        self.unalign_label2.configure(justify='left')
        self.unalign_label2.configure(text='''Unaligned FASTA folder''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_tooltip = ToolTip(self.unalign_label2, self.tooltip_font,
                                       '''unaligned''')
        self.fasta_folder_entry = ttk.Entry(self.Labelframe1)
        self.fasta_folder_entry.place(relx=0.314, rely=0.4, relheight=0.233
                                      , relwidth=0.489, bordermode='ignore')
        self.fasta_folder_entry.configure(textvariable=self.fasta_folder)
        self.fasta_folder_entry.configure(takefocus="")
        self.fasta_folder_entry.configure(cursor="fleur")
        self.fasta_folder_entry_tooltip = ToolTip(
            self.fasta_folder_entry, self.tooltip_font,
            'unaligned fasta files')
        self.open1_btn = tk.Button(self.Labelframe1)
        self.open1_btn.place(relx=0.82, rely=0.4, height=35, width=90
                             , bordermode='ignore')
        self.open1_btn.configure(activebackground="#ececec")
        self.open1_btn.configure(activeforeground="#000000")
        self.open1_btn.configure(background="#edf0f3")
        self.open1_btn.configure(command=open_file(self.fasta_folder_entry,
                                                   type_='folder'))
        self.open1_btn.configure(compound='left')
        self.open1_btn.configure(font="-family {TkDefaultFont} -size 10")
        self.open1_btn.configure(foreground="#000000")
        self.open1_btn.configure(highlightbackground="#edf0f3")
        self.open1_btn.configure(highlightcolor="black")
        self.open1_btn.configure(pady="0")
        self.open1_btn.configure(text='''Open''')

        self.align_label = tk.Label(self.Labelframe1)
        self.align_label.place(relx=0.03, rely=0.667, height=35, width=150
                               , bordermode='ignore')
        self.align_label.configure(activebackground="#f9f9f9")
        self.align_label.configure(activeforeground="black")
        self.align_label.configure(anchor='w')
        self.align_label.configure(background="#edf0f3")
        self.align_label.configure(compound='left')
        self.align_label.configure(font="-family {TkDefaultFont} -size 10")
        self.align_label.configure(foreground="#000000")
        self.align_label.configure(highlightbackground="#edf0f3")
        self.align_label.configure(highlightcolor="black")
        self.align_label.configure(justify='left')
        self.align_label.configure(text='''Aligned FASTA files''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_3_tooltip = ToolTip(self.align_label, self.tooltip_font,
                                         '''aligned''')

        self.aln_entry = ttk.Entry(self.Labelframe1)
        self.aln_entry.place(relx=0.314, rely=0.667, relheight=0.233
                             , relwidth=0.489, bordermode='ignore')
        self.aln_entry.configure(textvariable=self.aln)
        self.aln_entry.configure(takefocus="")
        self.aln_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.aln_entry_tooltip = ToolTip(self.aln_entry, self.tooltip_font,
                                         '''aligned fasta files''')
        self.open2_btn = tk.Button(self.Labelframe1)
        self.open2_btn.place(relx=0.82, rely=0.667, height=35, width=90
                             , bordermode='ignore')
        self.open2_btn.configure(activebackground="#ececec")
        self.open2_btn.configure(activeforeground="#000000")
        self.open2_btn.configure(background="#edf0f3")
        self.open2_btn.configure(command=open_file(self.aln_entry,
                                                   single=False))
        self.open2_btn.configure(compound='left')
        self.open2_btn.configure(font="-family {TkDefaultFont} -size 10")
        self.open2_btn.configure(foreground="#000000")
        self.open2_btn.configure(highlightbackground="#edf0f3")
        self.open2_btn.configure(highlightcolor="black")
        self.open2_btn.configure(pady="0")
        self.open2_btn.configure(text='''Open''')

        self.out_label = tk.Label(self.top)
        self.out_label.place(relx=0.117, rely=0.356, height=35, width=100)
        self.out_label.configure(activebackground="#f9f9f9")
        self.out_label.configure(activeforeground="black")
        self.out_label.configure(anchor='w')
        self.out_label.configure(background="#edf0f3")
        self.out_label.configure(compound='left')
        self.out_label.configure(font="-family {TkDefaultFont} -size 10")
        self.out_label.configure(foreground="#000000")
        self.out_label.configure(highlightbackground="#edf0f3")
        self.out_label.configure(highlightcolor="black")
        self.out_label.configure(justify='left')
        self.out_label.configure(text='''Output folder''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_3_1_tooltip = ToolTip(self.out_label, self.tooltip_font,
                                           '''output''')

        self.out_entry = ttk.Entry(self.top)
        self.out_entry.place(relx=0.317, rely=0.356, relheight=0.078
                             , relwidth=0.467)
        self.out_entry.configure(textvariable=self.out)
        self.out_entry.configure(takefocus="")
        self.out_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.out_entry_tooltip = \
            ToolTip(self.out_entry, self.tooltip_font,
                    '''unaligned fasta files''')

        self.open3_btn = tk.Button(self.top)
        self.open3_btn.place(relx=0.8, rely=0.356, height=35, width=90)
        self.open3_btn.configure(activebackground="#ececec")
        self.open3_btn.configure(activeforeground="#000000")
        self.open3_btn.configure(background="#edf0f3")
        self.open3_btn.configure(command=open_file(self.out_entry,
                                                   type_='folder'))
        self.open3_btn.configure(compound='left')
        self.open3_btn.configure(font="-family {TkDefaultFont} -size 10")
        self.open3_btn.configure(foreground="#000000")
        self.open3_btn.configure(highlightbackground="#edf0f3")
        self.open3_btn.configure(highlightcolor="black")
        self.open3_btn.configure(pady="0")
        self.open3_btn.configure(text='''Open''')

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.467, relheight=0.156
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Sliding window''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(cursor="fleur")
        self.Labelframe1.configure(highlightbackground="#d9d9d9")
        self.Labelframe1.configure(highlightcolor="black")

        self.Checkbutton1 = tk.Checkbutton(self.Labelframe1)
        self.Checkbutton1.place(relx=0.035, rely=0.429, relheight=0.3
                                , relwidth=0.262, bordermode='ignore')
        self.Checkbutton1.configure(activebackground="#ececec")
        self.Checkbutton1.configure(activeforeground="#000000")
        self.Checkbutton1.configure(anchor='w')
        self.Checkbutton1.configure(background="#edf0f3")
        self.Checkbutton1.configure(compound='left')
        self.Checkbutton1.configure(foreground="#000000")
        self.Checkbutton1.configure(highlightbackground="#edf0f3")
        self.Checkbutton1.configure(highlightcolor="black")
        self.Checkbutton1.configure(justify='left')
        self.Checkbutton1.configure(selectcolor="#edf0f3")
        self.Checkbutton1.configure(text='''Skip sliding window''')
        self.Checkbutton1.configure(variable=self.quick)

        self.window_size_label = tk.Label(self.Labelframe1)
        self.window_size_label.place(relx=0.314, rely=0.429, height=22, width=80
                                     , bordermode='ignore')
        self.window_size_label.configure(activebackground="#f9f9f9")
        self.window_size_label.configure(anchor='w')
        self.window_size_label.configure(background="#edf0f3")
        self.window_size_label.configure(compound='left')
        self.window_size_label.configure(foreground="#000000")
        self.window_size_label.configure(highlightbackground="#d9d9d9")
        self.window_size_label.configure(highlightcolor="black")
        self.window_size_label.configure(text='''Window size''')

        self.step_len_label = tk.Label(self.Labelframe1)
        self.step_len_label.place(relx=0.681, rely=0.429, height=22, width=80
                                  , bordermode='ignore')
        self.step_len_label.configure(activebackground="#f9f9f9")
        self.step_len_label.configure(anchor='w')
        self.step_len_label.configure(background="#edf0f3")
        self.step_len_label.configure(compound='left')
        self.step_len_label.configure(foreground="#000000")
        self.step_len_label.configure(highlightbackground="#d9d9d9")
        self.step_len_label.configure(highlightcolor="black")
        self.step_len_label.configure(text='''Step length''')

        self.size_entry = ttk.Entry(self.Labelframe1)
        self.size_entry.place(relx=0.489, rely=0.357, relheight=0.5
                              , relwidth=0.14, bordermode='ignore')
        self.size_entry.configure(textvariable=self.size)
        self.size_entry.configure(takefocus="")
        self.size_entry.configure(cursor="fleur")

        self.step_entry = ttk.Entry(self.Labelframe1)
        self.step_entry.place(relx=0.838, rely=0.357, relheight=0.5
                              , relwidth=0.14, bordermode='ignore')
        self.step_entry.configure(textvariable=self.step)
        self.step_entry.configure(takefocus="")
        self.step_entry.configure(cursor="fleur")
        self.size_entry.insert(0, '500')
        self.step_entry.insert(0, '50')

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.644, relheight=0.156
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Advance''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(highlightbackground="#d9d9d9")
        self.Labelframe1.configure(highlightcolor="black")

        self.Checkbutton1 = tk.Checkbutton(self.Labelframe1)
        self.Checkbutton1.place(relx=0.087, rely=0.429, relheight=0.3
                                , relwidth=0.314, bordermode='ignore')
        self.Checkbutton1.configure(activebackground="#ececec")
        self.Checkbutton1.configure(activeforeground="#000000")
        self.Checkbutton1.configure(anchor='w')
        self.Checkbutton1.configure(background="#edf0f3")
        self.Checkbutton1.configure(compound='left')
        self.Checkbutton1.configure(foreground="#000000")
        self.Checkbutton1.configure(highlightbackground="#edf0f3")
        self.Checkbutton1.configure(highlightcolor="black")
        self.Checkbutton1.configure(justify='left')
        self.Checkbutton1.configure(selectcolor="#edf0f3")
        self.Checkbutton1.configure(text='''Ignore gaps in alignment''')
        self.Checkbutton1.configure(variable=self.ig)

        self.Checkbutton1_2 = tk.Checkbutton(self.Labelframe1)
        self.Checkbutton1_2.place(relx=0.593, rely=0.429, relheight=0.3
                                  , relwidth=0.314, bordermode='ignore')
        self.Checkbutton1_2.configure(activebackground="#ececec")
        self.Checkbutton1_2.configure(activeforeground="#000000")
        self.Checkbutton1_2.configure(anchor='w')
        self.Checkbutton1_2.configure(background="#edf0f3")
        self.Checkbutton1_2.configure(compound='left')
        self.Checkbutton1_2.configure(foreground="#000000")
        self.Checkbutton1_2.configure(highlightbackground="#edf0f3")
        self.Checkbutton1_2.configure(highlightcolor="black")
        self.Checkbutton1_2.configure(justify='left')
        self.Checkbutton1_2.configure(selectcolor="#edf0f3")
        self.Checkbutton1_2.configure(text='''Ignore ambiguous bases''')
        self.Checkbutton1_2.configure(variable=self.iab)

        self.Button1_3 = tk.Button(self.top)
        self.Button1_3.place(relx=0.333, rely=0.867, height=40, width=189)
        self.Button1_3.configure(activebackground="#ececec")
        self.Button1_3.configure(activeforeground="#000000")
        self.Button1_3.configure(background="#edf0f3")
        self.Button1_3.configure(command=run_evaluate(self, self.top))
        self.Button1_3.configure(compound='left')
        self.Button1_3.configure(font="-family {TkDefaultFont} -size 12")
        self.Button1_3.configure(foreground="#000000")
        self.Button1_3.configure(highlightbackground="#edf0f3")
        self.Button1_3.configure(highlightcolor="black")
        self.Button1_3.configure(pady="0")
        self.Button1_3.configure(relief="raised")
        self.Button1_3.configure(text='''Run''')


class Primer:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#edf0f3'  # Closest X11 color: 'gray94'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9'  # X11 color: 'gray85'
        _ana1color = '#d9d9d9'  # X11 color: 'gray85'
        _ana2color = '#ececec'  # Closest X11 color: 'gray92'
        _tabfg1 = 'black'
        _tabfg2 = 'black'
        _tabbg1 = 'grey75'
        _tabbg2 = 'grey89'
        _bgmode = 'light'
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.', background=_bgcolor)
        self.style.configure('.', foreground=_fgcolor)
        self.style.configure('.', font="TkDefaultFont")
        self.style.map('.', background=
        [('selected', _compcolor), ('active', _ana2color)])

        top.geometry("600x500+284+458")
        move_to_center(top, 600, 600)
        top.title("Primer")
        top.configure(background="#edf0f3")
        top.configure(highlightbackground="#d9d9d9")
        top.configure(highlightcolor="black")

        self.top = top
        self.aln = tk.StringVar()
        self.aln_folder = tk.StringVar()
        self.out = tk.StringVar()
        self.coverage = tk.StringVar()
        self.mismatch = tk.StringVar()
        self.resolution = tk.StringVar()
        self.top_n = tk.StringVar()
        self.pmin = tk.StringVar()
        self.pmax = tk.StringVar()
        self.amin = tk.StringVar()
        self.amax = tk.StringVar()
        self.size = tk.StringVar()
        self.step = tk.StringVar()

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.02, relheight=0.2
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Input''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(highlightbackground="#d9d9d9")
        self.Labelframe1.configure(highlightcolor="black")

        self.TLabel1 = tk.Label(self.Labelframe1)
        self.TLabel1.place(relx=0.052, rely=0.2, height=35, width=130
                           , bordermode='ignore')
        self.TLabel1.configure(activebackground="#f9f9f9")
        self.TLabel1.configure(activeforeground="black")
        self.TLabel1.configure(anchor='w')
        self.TLabel1.configure(background="#edf0f3")
        self.TLabel1.configure(compound='left')
        self.TLabel1.configure(font="-family {TkDefaultFont} -size 10")
        self.TLabel1.configure(foreground="#000000")
        self.TLabel1.configure(highlightbackground="#edf0f3")
        self.TLabel1.configure(highlightcolor="black")
        self.TLabel1.configure(justify='left')
        self.TLabel1.configure(text='''Aligned FASTA files''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_tooltip = \
            ToolTip(self.TLabel1, self.tooltip_font, '''unaligned''')

        self.aln_entry = ttk.Entry(self.Labelframe1)
        self.aln_entry.place(relx=0.314, rely=0.2, relheight=0.35,
                             relwidth=0.489
                             , bordermode='ignore')
        self.aln_entry.configure(textvariable=self.aln)
        self.aln_entry.configure(takefocus="")
        self.aln_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.aln_entry_tooltip = \
            ToolTip(self.aln_entry, self.tooltip_font,
                    '''unaligned fasta files''')

        self.out_b = tk.Button(self.Labelframe1)
        self.out_b.place(relx=0.82, rely=0.2, height=35, width=90
                         , bordermode='ignore')
        self.out_b.configure(activebackground="#ececec")
        self.out_b.configure(activeforeground="#000000")
        self.out_b.configure(background="#edf0f3")
        self.out_b.configure(command=open_file(self.aln_entry, single=False))
        self.out_b.configure(compound='left')
        self.out_b.configure(font="-family {TkDefaultFont} -size 10")
        self.out_b.configure(foreground="#000000")
        self.out_b.configure(highlightbackground="#edf0f3")
        self.out_b.configure(highlightcolor="black")
        self.out_b.configure(pady="0")
        self.out_b.configure(relief="raised")
        self.out_b.configure(text='''Open''')

        self.TLabel1 = tk.Label(self.Labelframe1)
        self.TLabel1.place(relx=0.052, rely=0.6, height=35, width=160
                           , bordermode='ignore')
        self.TLabel1.configure(activebackground="#f9f9f9")
        self.TLabel1.configure(activeforeground="black")
        self.TLabel1.configure(anchor='w')
        self.TLabel1.configure(background="#edf0f3")
        self.TLabel1.configure(compound='left')
        self.TLabel1.configure(font="-family {TkDefaultFont} -size 10")
        self.TLabel1.configure(foreground="#000000")
        self.TLabel1.configure(highlightbackground="#edf0f3")
        self.TLabel1.configure(highlightcolor="black")
        self.TLabel1.configure(justify='left')
        self.TLabel1.configure(text='''Aligned FASTA folder''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_tooltip = \
            ToolTip(self.TLabel1, self.tooltip_font, '''unaligned''')

        self.aln_folder_entry = ttk.Entry(self.Labelframe1)
        self.aln_folder_entry.place(relx=0.314, rely=0.6, relheight=0.35
                                    , relwidth=0.489, bordermode='ignore')
        self.aln_folder_entry.configure(textvariable=self.aln_folder)
        self.aln_folder_entry.configure(takefocus="")
        self.aln_folder_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.aln_folder_entry_tooltip = \
            ToolTip(self.aln_folder_entry, self.tooltip_font,
                    '''unaligned fasta files''')

        self.folder_b = tk.Button(self.Labelframe1)
        self.folder_b.place(relx=0.82, rely=0.6, height=35, width=90
                            , bordermode='ignore')
        self.folder_b.configure(activebackground="#ececec")
        self.folder_b.configure(activeforeground="#000000")
        self.folder_b.configure(background="#edf0f3")
        self.folder_b.configure(command=open_file(self.aln_folder_entry,
                                                  type_='folder'))
        self.folder_b.configure(compound='left')
        self.folder_b.configure(font="-family {TkDefaultFont} -size 10")
        self.folder_b.configure(foreground="#000000")
        self.folder_b.configure(highlightbackground="#edf0f3")
        self.folder_b.configure(highlightcolor="black")
        self.folder_b.configure(pady="0")
        self.folder_b.configure(relief="raised")
        self.folder_b.configure(text='''Open''')

        self.out_label = tk.Label(self.top)
        self.out_label.place(relx=0.07, rely=0.24, height=35, width=150)
        self.out_label.configure(activebackground="#f9f9f9")
        self.out_label.configure(activeforeground="black")
        self.out_label.configure(anchor='w')
        self.out_label.configure(background="#edf0f3")
        self.out_label.configure(compound='left')
        self.out_label.configure(font="-family {TkDefaultFont} -size 10")
        self.out_label.configure(foreground="#000000")
        self.out_label.configure(highlightbackground="#edf0f3")
        self.out_label.configure(highlightcolor="black")
        self.out_label.configure(justify='left')
        self.out_label.configure(text='''Output folder''')
        self.tooltip_font = "TkDefaultFont"
        self.TLabel1_3_1_tooltip = \
            ToolTip(self.out_label, self.tooltip_font, '''Output''')

        self.out_entry = ttk.Entry(self.top)
        self.out_entry.place(relx=0.317, rely=0.24, relheight=0.07
                             , relwidth=0.467)
        self.out_entry.configure(textvariable=self.out)
        self.out_entry.configure(takefocus="")
        self.out_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.out_entry_tooltip = \
            ToolTip(self.out_entry, self.tooltip_font,
                    '''unaligned fasta files''')

        self.out_b = tk.Button(self.top)
        self.out_b.place(relx=0.8, rely=0.24, height=35, width=90)
        self.out_b.configure(activebackground="#ececec")
        self.out_b.configure(activeforeground="#000000")
        self.out_b.configure(background="#edf0f3")
        self.out_b.configure(command=open_file(self.out_entry, type_='folder'))
        self.out_b.configure(compound='left')
        self.out_b.configure(font="-family {TkDefaultFont} -size 10")
        self.out_b.configure(foreground="#000000")
        self.out_b.configure(highlightbackground="#edf0f3")
        self.out_b.configure(highlightcolor="black")
        self.out_b.configure(pady="0")
        self.out_b.configure(relief="raised")
        self.out_b.configure(text='''Open''')

        self.Labelframe1 = tk.LabelFrame(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.34, relheight=0.38
                               , relwidth=0.955)
        self.Labelframe1.configure(relief='groove')
        self.Labelframe1.configure(foreground="#000000")
        self.Labelframe1.configure(text='''Advance''')
        self.Labelframe1.configure(background="#edf0f3")
        self.Labelframe1.configure(highlightbackground="#edf0f3")
        self.Labelframe1.configure(highlightcolor="black")

        self.coverage_label = tk.Label(self.Labelframe1)
        self.coverage_label.place(relx=0.07, rely=0.158, height=35
                                  , width=60, bordermode='ignore')
        self.coverage_label.configure(activebackground="#f9f9f9")
        self.coverage_label.configure(activeforeground="black")
        self.coverage_label.configure(anchor='w')
        self.coverage_label.configure(background="#edf0f3")
        self.coverage_label.configure(compound='left')
        self.coverage_label.configure(font="-family {TkDefaultFont} -size 10")
        self.coverage_label.configure(foreground="#000000")
        self.coverage_label.configure(highlightbackground="#edf0f3")
        self.coverage_label.configure(highlightcolor="black")
        self.coverage_label.configure(justify='left')
        self.coverage_label.configure(text='''Coverage''')

        self.coverage_entry = ttk.Entry(self.Labelframe1)
        self.coverage_entry.place(relx=0.297, rely=0.142, relheight=0.184
                                  , relwidth=0.192, bordermode='ignore')
        self.coverage_entry.configure(textvariable=self.coverage)
        self.coverage_entry.configure(takefocus="")
        self.coverage_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.coverage_entry_tooltip = \
            ToolTip(self.coverage_entry, self.tooltip_font,
                    '''minimal coverage of primer on alignment''')
        self.coverage_entry.insert(0, '0.5')

        self.mismatch_label = tk.Label(self.Labelframe1)
        self.mismatch_label.place(relx=0.07, rely=0.368, height=35
                                  , width=120, bordermode='ignore')
        self.mismatch_label.configure(activebackground="#f9f9f9")
        self.mismatch_label.configure(activeforeground="black")
        self.mismatch_label.configure(anchor='w')
        self.mismatch_label.configure(background="#edf0f3")
        self.mismatch_label.configure(compound='left')
        self.mismatch_label.configure(font="-family {TkDefaultFont} -size 10")
        self.mismatch_label.configure(foreground="#000000")
        self.mismatch_label.configure(highlightbackground="#edf0f3")
        self.mismatch_label.configure(highlightcolor="black")
        self.mismatch_label.configure(justify='left')
        self.mismatch_label.configure(text='''Mismatch''')

        self.res_label = tk.Label(self.Labelframe1)
        self.res_label.place(relx=0.07, rely=0.579, height=35
                             , width=130, bordermode='ignore')
        self.res_label.configure(activebackground="#f9f9f9")
        self.res_label.configure(activeforeground="black")
        self.res_label.configure(anchor='w')
        self.res_label.configure(background="#edf0f3")
        self.res_label.configure(compound='left')
        self.res_label.configure(font="-family {TkDefaultFont} -size 10")
        self.res_label.configure(foreground="#000000")
        self.res_label.configure(highlightbackground="#edf0f3")
        self.res_label.configure(highlightcolor="black")
        self.res_label.configure(justify='left')
        self.res_label.configure(text='''Resolution''')

        self.mismatch_entry = ttk.Entry(self.Labelframe1)
        self.mismatch_entry.place(relx=0.297, rely=0.368, relheight=0.184
                                  , relwidth=0.192, bordermode='ignore')
        self.mismatch_entry.configure(textvariable=self.mismatch)
        self.mismatch_entry.configure(takefocus="")
        self.mismatch_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.mismatch_entry_tooltip = \
            ToolTip(self.mismatch_entry, self.tooltip_font,
                    '''maximum mismatch bases in primer''')

        self.resolution_entry = ttk.Entry(self.Labelframe1)
        self.resolution_entry.place(relx=0.297, rely=0.579, relheight=0.184
                                    , relwidth=0.192, bordermode='ignore')
        self.resolution_entry.configure(textvariable=self.resolution)
        self.resolution_entry.configure(takefocus="")
        self.resolution_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.resolution_entry_tooltip = \
            ToolTip(self.resolution_entry, self.tooltip_font,
                    '''minimal resolution of amplified fragment''')
        self.mismatch_entry.insert(0, '4')
        self.resolution_entry.insert(0, '0.3')

        self.TSeparator2 = ttk.Separator(self.Labelframe1)
        self.TSeparator2.place(relx=0.524, rely=0.137, relheight=0.779
                               , bordermode='ignore')
        self.TSeparator2.configure(orient="vertical")

        self.topn_label = tk.Label(self.Labelframe1)
        self.topn_label.place(relx=0.07, rely=0.789, height=35
                              , width=130, bordermode='ignore')
        self.topn_label.configure(activebackground="#f9f9f9")
        self.topn_label.configure(activeforeground="black")
        self.topn_label.configure(anchor='w')
        self.topn_label.configure(background="#edf0f3")
        self.topn_label.configure(compound='left')
        self.topn_label.configure(font="-family {TkDefaultFont} -size 10")
        self.topn_label.configure(foreground="#000000")
        self.topn_label.configure(highlightbackground="#edf0f3")
        self.topn_label.configure(highlightcolor="black")
        self.topn_label.configure(justify='left')
        self.topn_label.configure(text='''Top n''')

        self.top_n_entry = ttk.Entry(self.Labelframe1)
        self.top_n_entry.place(relx=0.297, rely=0.789, relheight=0.184
                               , relwidth=0.192, bordermode='ignore')
        self.top_n_entry.configure(textvariable=self.top_n)
        self.top_n_entry.configure(takefocus="")
        self.top_n_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.top_n_entry_tooltip = \
            ToolTip(self.top_n_entry, self.tooltip_font,
                    '''Only keep top best primers''')
        self.top_n_entry.insert(0, '1')

        self.primer_len_label = tk.Label(self.Labelframe1)
        self.primer_len_label.place(relx=0.541, rely=0.158, height=35
                                    , width=100, bordermode='ignore')
        self.primer_len_label.configure(activebackground="#f9f9f9")
        self.primer_len_label.configure(activeforeground="black")
        self.primer_len_label.configure(anchor='w')
        self.primer_len_label.configure(background="#edf0f3")
        self.primer_len_label.configure(compound='left')
        self.primer_len_label.configure(font="-family {TkDefaultFont} -size 10")
        self.primer_len_label.configure(foreground="#000000")
        self.primer_len_label.configure(highlightbackground="#edf0f3")
        self.primer_len_label.configure(highlightcolor="black")
        self.primer_len_label.configure(justify='left')
        self.primer_len_label.configure(text='''Primer length''')

        self.pmin_entry = ttk.Entry(self.Labelframe1)
        self.pmin_entry.place(relx=0.716, rely=0.158, relheight=0.184
                              , relwidth=0.087, bordermode='ignore')
        self.pmin_entry.configure(textvariable=self.pmin)
        self.pmin_entry.configure(takefocus="")
        self.pmin_entry.configure(cursor="fleur")

        self.Label1 = tk.Label(self.Labelframe1)
        self.Label1.place(relx=0.803, rely=0.158, height=35, width=36
                          , bordermode='ignore')
        self.Label1.configure(activebackground="#f9f9f9")
        self.Label1.configure(activeforeground="SystemButtonText")
        self.Label1.configure(anchor='w')
        self.Label1.configure(background="#edf0f3")
        self.Label1.configure(compound='left')
        self.Label1.configure(font="-family {TkDefaultFont} -size 13")
        self.Label1.configure(foreground="#000000")
        self.Label1.configure(highlightbackground="#edf0f3")
        self.Label1.configure(highlightcolor="black")
        self.Label1.configure(text='''to''')

        self.pmax_entry = ttk.Entry(self.Labelframe1)
        self.pmax_entry.place(relx=0.855, rely=0.158, relheight=0.184
                              , relwidth=0.122, bordermode='ignore')
        self.pmax_entry.configure(textvariable=self.pmax)
        self.pmax_entry.configure(takefocus="")
        self.pmax_entry.configure(cursor="fleur")
        self.pmin_entry.insert(0, '20')
        self.pmax_entry.insert(0, '30')

        self.amp_len_label = tk.Label(self.Labelframe1)
        self.amp_len_label.place(relx=0.541, rely=0.368, height=35
                                 , width=120, bordermode='ignore')
        self.amp_len_label.configure(activebackground="#f9f9f9")
        self.amp_len_label.configure(activeforeground="black")
        self.amp_len_label.configure(anchor='w')
        self.amp_len_label.configure(background="#edf0f3")
        self.amp_len_label.configure(compound='left')
        self.amp_len_label.configure(font="-family {TkDefaultFont} -size 10")
        self.amp_len_label.configure(foreground="#000000")
        self.amp_len_label.configure(highlightbackground="#edf0f3")
        self.amp_len_label.configure(highlightcolor="black")
        self.amp_len_label.configure(justify='left')
        self.amp_len_label.configure(text='''Amplicon size''')

        self.amin_entry = ttk.Entry(self.Labelframe1)
        self.amin_entry.place(relx=0.716, rely=0.368, relheight=0.184
                              , relwidth=0.087, bordermode='ignore')
        self.amin_entry.configure(textvariable=self.amin)
        self.amin_entry.configure(takefocus="")
        self.amin_entry.configure(cursor="fleur")
        self.tooltip_font = "TkDefaultFont"
        self.amin_entry_tooltip = \
            ToolTip(self.amin_entry, self.tooltip_font,
                    '''including primer length''')

        self.Label1 = tk.Label(self.Labelframe1)
        self.Label1.place(relx=0.803, rely=0.368, height=35, width=36
                          , bordermode='ignore')
        self.Label1.configure(activebackground="#f9f9f9")
        self.Label1.configure(activeforeground="SystemButtonText")
        self.Label1.configure(anchor='w')
        self.Label1.configure(background="#edf0f3")
        self.Label1.configure(compound='left')
        self.Label1.configure(font="-family {TkDefaultFont} -size 13")
        self.Label1.configure(foreground="#000000")
        self.Label1.configure(highlightbackground="#edf0f3")
        self.Label1.configure(highlightcolor="black")
        self.Label1.configure(text='''to''')

        self.amax_entry = ttk.Entry(self.Labelframe1)
        self.amax_entry.place(relx=0.855, rely=0.368, relheight=0.184
                              , relwidth=0.122, bordermode='ignore')
        self.amax_entry.configure(textvariable=self.amax)
        self.amax_entry.configure(takefocus="")
        self.amax_entry.configure(cursor="fleur")
        self.amin_entry.insert(0, '300')
        self.amax_entry.insert(0, '800')

        self.sliding_window_label = tk.Label(self.Labelframe1)
        self.sliding_window_label.place(relx=0.541, rely=0.579, height=35,
                                        width=130
                                        , bordermode='ignore')
        self.sliding_window_label.configure(activebackground="#f9f9f9")
        self.sliding_window_label.configure(anchor='w')
        self.sliding_window_label.configure(background="#edf0f3")
        self.sliding_window_label.configure(compound='left')
        self.sliding_window_label.configure(foreground="#000000")
        self.sliding_window_label.configure(highlightbackground="#d9d9d9")
        self.sliding_window_label.configure(highlightcolor="black")
        self.sliding_window_label.configure(text='''Sliding window size''')

        self.sliding_window2_label = tk.Label(self.Labelframe1)
        self.sliding_window2_label.place(relx=0.541, rely=0.789, height=35,
                                         width=130
                                         , bordermode='ignore')
        self.sliding_window2_label.configure(activebackground="#f9f9f9")
        self.sliding_window2_label.configure(anchor='w')
        self.sliding_window2_label.configure(background="#edf0f3")
        self.sliding_window2_label.configure(compound='left')
        self.sliding_window2_label.configure(foreground="#000000")
        self.sliding_window2_label.configure(highlightbackground="#d9d9d9")
        self.sliding_window2_label.configure(highlightcolor="black")
        self.sliding_window2_label.configure(text='''Sliding window step''')

        self.size_entry = ttk.Entry(self.Labelframe1)
        self.size_entry.place(relx=0.803, rely=0.579, relheight=0.184
                              , relwidth=0.14, bordermode='ignore')
        self.size_entry.configure(textvariable=self.size)
        self.size_entry.configure(takefocus="")
        self.size_entry.configure(cursor="fleur")

        self.step_entry = ttk.Entry(self.Labelframe1)
        self.step_entry.place(relx=0.803, rely=0.789, relheight=0.184
                              , relwidth=0.14, bordermode='ignore')
        self.step_entry.configure(textvariable=self.step)
        self.step_entry.configure(takefocus="")
        self.step_entry.configure(cursor="fleur")
        self.size_entry.insert(0, '500')
        self.step_entry.insert(0, '50')

        self.run_b = tk.Button(self.top)
        self.run_b.place(relx=0.367, rely=0.82, height=40, width=180)
        self.run_b.configure(activebackground="#ececec")
        self.run_b.configure(activeforeground="#000000")
        self.run_b.configure(background="#edf0f3")
        self.run_b.configure(command=run_primer(self, self.top))
        self.run_b.configure(compound='left')
        self.run_b.configure(font="-family {TkDefaultFont} -size 14")
        self.run_b.configure(foreground="#000000")
        self.run_b.configure(highlightbackground="#edf0f3")
        self.run_b.configure(highlightcolor="black")
        self.run_b.configure(pady="0")
        self.run_b.configure(relief="raised")
        self.run_b.configure(text='''Run''')


class LOG:
    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#edf0f3'  # Closest X11 color: 'gray94'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9'  # X11 color: 'gray85'
        _ana1color = '#d9d9d9'  # X11 color: 'gray85'
        _ana2color = '#ececec'  # Closest X11 color: 'gray92'
        _tabfg1 = 'black'
        _tabfg2 = 'black'
        _tabbg1 = 'grey75'
        _tabbg2 = 'grey89'
        _bgmode = 'light'
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.', background=_bgcolor)
        self.style.configure('.', foreground=_fgcolor)
        self.style.map('.', background=
        [('selected', _compcolor), ('active', _ana2color)])

        top.geometry("600x400+2532+335")
        top.minsize(72, 15)
        top.maxsize(3648, 1052)
        top.resizable(1, 1)
        top.title("Running...")
        top.configure(background="#edf0f3")
        top.configure(highlightbackground="#d9d9d9")
        top.configure(highlightcolor="black")

        self.top = top

        self.Scrolledtext1 = ScrolledText(self.top)
        self.Scrolledtext1.place(relx=0.0, rely=0.0, relheight=1.0,
                                 relwidth=1.0)

        self.Scrolledtext1.configure(background="white")
        self.Scrolledtext1.configure(font="TkTextFont")
        self.Scrolledtext1.configure(foreground="black")
        self.Scrolledtext1.configure(highlightbackground="#edf0f3")
        self.Scrolledtext1.configure(highlightcolor="black")
        self.Scrolledtext1.configure(insertbackground="black")
        self.Scrolledtext1.configure(insertborderwidth="3")
        self.Scrolledtext1.configure(selectbackground="#c4c4c4")
        self.Scrolledtext1.configure(selectforeground="black")
        self.Scrolledtext1.configure(wrap="none")


class ToolTip(tk.Toplevel):
    """ Provides a ToolTip widget for Tkinter. """

    def __init__(self, wdgt, tooltip_font, msg=None, msgFunc=None,
                 delay=0.5, follow=True):
        self.wdgt = wdgt
        self.parent = self.wdgt.master
        tk.Toplevel.__init__(self, self.parent, bg='black', padx=1, pady=1)
        self.withdraw()
        self.overrideredirect(True)
        self.msgVar = tk.StringVar()
        if msg is None:
            self.msgVar.set('No message provided')
        else:
            self.msgVar.set(msg)
        self.msgFunc = msgFunc
        self.delay = delay
        self.follow = follow
        self.visible = 0
        self.lastMotion = 0
        tk.Message(self, textvariable=self.msgVar, bg='#FFFFDD',
                   font=tooltip_font,
                   aspect=1000).grid()
        self.wdgt.bind('<Enter>', self.spawn, '+')
        self.wdgt.bind('<Leave>', self.hide, '+')
        self.wdgt.bind('<Motion>', self.move, '+')

    def spawn(self, event=None):
        self.visible = 1
        self.after(int(self.delay * 1000), self.show)

    def show(self):
        if self.visible == 1 and time() - self.lastMotion > self.delay:
            self.visible = 2
        if self.visible == 2:
            self.deiconify()

    def move(self, event):
        self.lastMotion = time()
        if self.follow is False:
            self.withdraw()
            self.visible = 1
        self.geometry('+%i+%i' % (event.x_root + 20, event.y_root - 10))
        try:
            self.msgVar.set(self.msgFunc())
        except:
            pass
        self.after(int(self.delay * 1000), self.show)

    def hide(self, event=None):
        self.visible = 0
        self.withdraw()

    def update(self, msg):
        self.msgVar.set(msg)


class AutoScroll(object):
    '''Configure the scrollbars for a widget.'''

    def __init__(self, master):
        #  Rozen. Added the try-except clauses so that this class
        #  could be used for scrolled entry widget for which vertical
        #  scrolling is not supported. 5/7/14.
        try:
            vsb = ttk.Scrollbar(master, orient='vertical', command=self.yview)
        except:
            pass
        hsb = ttk.Scrollbar(master, orient='horizontal', command=self.xview)
        try:
            self.configure(yscrollcommand=self._autoscroll(vsb))
        except:
            pass
        self.configure(xscrollcommand=self._autoscroll(hsb))
        self.grid(column=0, row=0, sticky='nsew')
        try:
            vsb.grid(column=1, row=0, sticky='ns')
        except:
            pass
        hsb.grid(column=0, row=1, sticky='ew')
        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)
        # Copy geometry methods of master  (taken from ScrolledText.py)
        methods = tk.Pack.__dict__.keys() | tk.Grid.__dict__.keys() \
                  | tk.Place.__dict__.keys()
        for meth in methods:
            if meth[0] != '_' and meth not in ('config', 'configure'):
                setattr(self, meth, getattr(master, meth))

    @staticmethod
    def _autoscroll(sbar):
        '''Hide and show scrollbar as needed.'''

        def wrapped(first, last):
            first, last = float(first), float(last)
            if first <= 0 and last >= 1:
                sbar.grid_remove()
            else:
                sbar.grid()
            sbar.set(first, last)

        return wrapped

    def __str__(self):
        return str(self.master)


def _create_container(func):
    '''Creates a ttk Frame with a given master, and use this new frame to
    place the scrollbars and the widget.'''

    def wrapped(cls, master, **kw):
        container = ttk.Frame(master)
        container.bind('<Enter>', lambda e: _bound_to_mousewheel(e, container))
        container.bind('<Leave>',
                       lambda e: _unbound_to_mousewheel(e, container))
        return func(cls, container, **kw)

    return wrapped


class ScrolledText(AutoScroll, tk.Text):
    '''A standard Tkinter Text widget with scrollbars that will
    automatically show/hide as needed.'''

    @_create_container
    def __init__(self, master, **kw):
        tk.Text.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)


def _bound_to_mousewheel(event, widget):
    child = widget.winfo_children()[0]
    if platform.system() == 'Windows' or platform.system() == 'Darwin':
        child.bind_all('<MouseWheel>', lambda e: _on_mousewheel(e, child))
        child.bind_all('<Shift-MouseWheel>', lambda e: _on_shiftmouse(e, child))
    else:
        child.bind_all('<Button-4>', lambda e: _on_mousewheel(e, child))
        child.bind_all('<Button-5>', lambda e: _on_mousewheel(e, child))
        child.bind_all('<Shift-Button-4>', lambda e: _on_shiftmouse(e, child))
        child.bind_all('<Shift-Button-5>', lambda e: _on_shiftmouse(e, child))


def _unbound_to_mousewheel(event, widget):
    if platform.system() == 'Windows' or platform.system() == 'Darwin':
        widget.unbind_all('<MouseWheel>')
        widget.unbind_all('<Shift-MouseWheel>')
    else:
        widget.unbind_all('<Button-4>')
        widget.unbind_all('<Button-5>')
        widget.unbind_all('<Shift-Button-4>')
        widget.unbind_all('<Shift-Button-5>')


def _on_mousewheel(event, widget):
    if platform.system() == 'Windows':
        widget.yview_scroll(-1 * int(event.delta / 120), 'units')
    elif platform.system() == 'Darwin':
        widget.yview_scroll(-1 * int(event.delta), 'units')
    else:
        if event.num == 4:
            widget.yview_scroll(-1, 'units')
        elif event.num == 5:
            widget.yview_scroll(1, 'units')


def _on_shiftmouse(event, widget):
    if platform.system() == 'Windows':
        widget.xview_scroll(-1 * int(event.delta / 120), 'units')
    elif platform.system() == 'Darwin':
        widget.xview_scroll(-1 * int(event.delta), 'units')
    else:
        if event.num == 4:
            widget.xview_scroll(-1, 'units')
        elif event.num == 5:
            widget.xview_scroll(1, 'units')


def ui_main():
    global root
    root = tk.Tk()
    root.protocol('WM_DELETE_WINDOW', root.destroy)
    # Creates a toplevel widget.
    global _top1, _w1
    _top1 = root
    _w1 = Root(_top1)
    root.mainloop()


def ui_gb2fasta():
    global _top2, _w2
    _top2 = tk.Toplevel(root)
    root.iconify()
    _top2.protocol('WM_DELETE_WINDOW', after_close(_top2))
    _w2 = GB2Fasta(_top2)


def ui_evaluate():
    global _top3, _w3
    _top3 = tk.Toplevel(root)
    root.iconify()
    _top3.protocol('WM_DELETE_WINDOW', after_close(_top3))
    _w3 = Evaluate(_top3)


def ui_primer():
    # Creates a toplevel widget.
    global _top4, _w4
    _top4 = tk.Toplevel(root)
    root.iconify()
    _top4.protocol('WM_DELETE_WINDOW', after_close(_top4))
    _w4 = Primer(_top4)


def get_arg_str(value: tk.Variable, name: str, arg_str: str,
                is_bool=False) -> str:
    value_str = ''
    if value.get():
        if is_bool:
            value_str = f'{name} '
        else:
            value_str = f'{name} {value.get()} '
    arg_str += value_str
    print(name, value_str)
    return arg_str


def run_gb2fasta(w: tk.Frame, t: tk.Toplevel):
    # todo: test options and functions
    def f():
        nonlocal w
        arg_str = ''
        arg_str = get_arg_str(w.gb, '-gb', arg_str)
        arg_str = get_arg_str(w.gene, '-gene', arg_str)
        arg_str = get_arg_str(w.molecular, '-molecular', arg_str)
        arg_str = get_arg_str(w.group, '-group', arg_str)
        arg_str = get_arg_str(w.og, '-og', arg_str)
        arg_str = get_arg_str(w.refseq, '-refseq', arg_str)
        arg_str = get_arg_str(w.count, '-count', arg_str)
        arg_str = get_arg_str(w.min_len, '-min_len', arg_str)
        arg_str = get_arg_str(w.max_len, '-max_len', arg_str)
        arg_str = get_arg_str(w.date_start, '-date_start', arg_str)
        arg_str = get_arg_str(w.date_end, '-date_end', arg_str)
        arg_str = get_arg_str(w.exclude, '-exclude', arg_str)
        arg_str = get_arg_str(w.query, '-query', arg_str)
        arg_str = get_arg_str(w.taxon, '-taxon', arg_str)
        arg_str = get_arg_str(w.out, '-out', arg_str)
        arg_str = get_arg_str(w.expand, '-expand', arg_str)
        arg_str = get_arg_str(w.max_name_len, '-max_name_len', arg_str)
        arg_str = get_arg_str(w.max_gene_len, '-max_gene_len', arg_str)
        arg_str = get_arg_str(w.unique, '-unique', arg_str)
        arg_str = get_arg_str(w.allow_repeat, '-allow_repeat', arg_str, is_bool=True)
        arg_str = get_arg_str(w.allow_invert_repeat, '-allow_invert_repeat', arg_str,
                    is_bool=True)
        arg_str = get_arg_str(w.allow_mosaic_repeat, '-allow_mosaic_repeat', arg_str,
                    is_bool=True)
        arg_str = get_arg_str(w.no_divide, '-no_divide', arg_str, is_bool=True)
        arg_str = get_arg_str(w.rename, '-rename', arg_str, is_bool=True)
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s*0.618)}+{w//3}+{h//3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        print(arg_str)
        r = threading.Thread(target=thread_wrap,
                             args=(gb2fasta_main, arg_str, run),
                             daemon=True)
        r.start()

    return f


def run_evaluate(w: tk.Frame, t: tk.Toplevel):
    # todo: test options and functions
    def f():
        nonlocal w
        arg_str = ''
        arg_str = get_arg_str(w.fasta, '-fasta', arg_str)
        arg_str = get_arg_str(w.fasta_folder, '-fasta_folder', arg_str)
        arg_str = get_arg_str(w.aln, '-aln', arg_str)
        arg_str = get_arg_str(w.out, '-out', arg_str)
        arg_str = get_arg_str(w.size, '-size', arg_str)
        arg_str = get_arg_str(w.step, '-step', arg_str)
        arg_str = get_arg_str(w.quick, '-quick', arg_str, is_bool=True)
        arg_str = get_arg_str(w.ig, '-ig', arg_str, is_bool=True)
        arg_str = get_arg_str(w.iab, '-iab', arg_str, is_bool=True)
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s*0.618)}+{w//3}+{h//3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        print(arg_str)
        r = threading.Thread(target=thread_wrap,
                             args=(evaluate_main, arg_str, run),
                             daemon=True)
        r.start()

    return f


def run_primer(w: tk.Frame, t: tk.Toplevel):
    # todo: test options and functions
    def f():
        nonlocal w
        arg_str = ''
        arg_str = get_arg_str(w.aln, '-aln', arg_str)
        arg_str = get_arg_str(w.aln_folder, '-aln_folder', arg_str)
        arg_str = get_arg_str(w.out, '-out', arg_str)
        arg_str = get_arg_str(w.coverage, '-coverage', arg_str)
        arg_str = get_arg_str(w.mismatch, '-mismatch', arg_str)
        arg_str = get_arg_str(w.resolution, '-resolution', arg_str)
        arg_str = get_arg_str(w.top_n, '-top_n', arg_str)
        arg_str = get_arg_str(w.pmin, '-pmin', arg_str)
        arg_str = get_arg_str(w.pmax, '-pmax', arg_str)
        arg_str = get_arg_str(w.amin, '-amin', arg_str)
        arg_str = get_arg_str(w.amax, '-amax', arg_str)
        arg_str = get_arg_str(w.size, '-size', arg_str)
        arg_str = get_arg_str(w.step, '-step', arg_str)
        arg_str += '-primer'
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s*0.618)}+{w//3}+{h//3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        print(arg_str)
        r = threading.Thread(target=thread_wrap,
                             args=(primer_main, arg_str, run),
                             daemon=True)
        r.start()

    return f


def run_help():
    url = 'https://github.com/wpwupingwp/barcodefinder'
    webbrowser.open(url, new=2)


if __name__ == '__main__':
    ui_main()