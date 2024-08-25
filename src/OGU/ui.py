#!/usr/bin/env python
import logging
import os
import platform
import queue
import sys
import threading
import tkinter as tk
import tkinter.ttk as ttk
import webbrowser
from importlib import resources
from logging import handlers
from time import time
from tkinter import filedialog, messagebox, scrolledtext

from OGU.evaluate import evaluate_main
from OGU.gb2fasta import gb2fasta_main
from OGU.global_vars import DATEFMT, FMT, log, name
from OGU.primer import primer_main
from OGU.utils import font_family, get_all_third_party


def my_labelframe(parent: tk.LabelFrame) -> tk.LabelFrame:
    frame = tk.LabelFrame(parent)
    frame.configure(background="#edf0f3")
    frame.configure(cursor="fleur")
    frame.configure(font=f'-family {font_family} -size 14')
    frame.configure(foreground="#000000")
    frame.configure(highlightbackground="#edf0f3")
    frame.configure(highlightcolor="black")
    frame.configure(relief='groove')
    return frame


def my_label(frame: tk.LabelFrame, fontsize=11) -> tk.Label:
    label = tk.Label(frame)
    label.configure(activebackground="#f9f9f9")
    label.configure(activeforeground="black")
    label.configure(anchor='w')
    label.configure(background="#edf0f3")
    label.configure(compound='left')
    label.configure(foreground="#000000")
    label.configure(highlightbackground="#edf0f3")
    label.configure(highlightcolor="black")
    label.configure(justify='left')
    label.configure(font=f'-family {font_family} -size {fontsize}')
    return label


def my_button(frame: tk.LabelFrame, fontsize=12) -> tk.Button:
    button = tk.Button(frame)
    button.configure(compound='left')
    button.configure(activebackground="#ececec")
    button.configure(activeforeground="#000000")
    button.configure(background="#edf0f3")
    button.configure(font=f"-family {font_family} -size {fontsize}")
    button.configure(foreground="#000000")
    button.configure(highlightbackground="#edf0f3")
    button.configure(highlightcolor="black")
    button.configure(pady="0")
    button.configure(relief="raised")
    button.configure(borderwidth=1)
    return button


def my_checkbutton(frame: tk.LabelFrame, fontsize=11) -> tk.Checkbutton:
    check_b = tk.Checkbutton(frame)
    check_b.configure(activebackground="#ececec")
    check_b.configure(activeforeground="#000000")
    check_b.configure(anchor='w')
    check_b.configure(background="#edf0f3")
    check_b.configure(compound='left')
    check_b.configure(foreground="#000000")
    check_b.configure(highlightbackground="#edf0f3")
    check_b.configure(highlightcolor="black")
    check_b.configure(justify='left')
    check_b.configure(selectcolor="#ffffff")
    check_b.configure(font=f'-family {font_family} -size {fontsize}')
    return check_b


class EntryWithPlaceholder(tk.Entry):
    def __init__(self, parent, placeholder='eg. ', color='grey'):
        super().__init__(parent, bg='white')
        self.placeholder = placeholder
        self.placeholder_color = color
        self.default_fg = self['fg']
        self.bind('<FocusOut>', self.focus_out)
        self.bind('<FocusIn>', self.focus_in)
        # tkinter handle a <FocusOut> immediately after entry created and
        # before the self.add_placeholder() line is executed
        self.after_idle(self.add_placeholder)

    def get(self):
        value = super().get()
        if value == self.placeholder:
            return ''
        else:
            return value

    def add_placeholder(self):
        self.config(fg=self.placeholder_color)
        self.insert(0, self.placeholder)

    def focus_in(self, *args):
        # remove placeholder when click
        if not self.get():
            self.delete('0', 'end')
            self.config(fg=self.default_fg)

    def focus_out(self, *args):
        # re-add placeholder if empty
        if not self.get():
            self.add_placeholder()


def my_entry(frame: tk.LabelFrame, placeholder='', color='grey',
             fontsize=11) -> tk.Entry:
    entry = EntryWithPlaceholder(frame, placeholder=placeholder, color=color)
    entry.configure(font=f'-family {font_family} -size {fontsize}')
    entry.configure(takefocus="")
    entry.configure(cursor="fleur")
    return entry


def my_combobox(frame: tk.LabelFrame, fontsize=11) -> ttk.Combobox:
    style = 'combostyle'
    combo_style = ttk.Style()
    if style not in ttk.Style().theme_names():
        combo_style.theme_create(style, parent='alt', settings={
            'TCombobox': {'configure': {
                'selectbackground': 'gray',
                'fieldbackground': 'white',
                'background': 'white'}}})
    combo_style.theme_use(style)
    combobox = ttk.Combobox(frame, font=(font_family, fontsize))
    combobox.configure(state='readonly')
    combobox.configure(takefocus="")
    return combobox


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


def thread_wrap(function, arg_str, window, no_arg=False):
    """
    Wrap for callback.
    Args:
        function(callable): function to call
        arg_str(str): string for function's argparse
        window(Toplevel): window to hide
        no_arg: function do not need args
    """
    try:
        if not no_arg:
            result = function(arg_str)
        else:
            result = function()
    except Exception as e:
        log.exception(str(e))
        log.exception('Abort.')
        messagebox.showinfo(message='Abort.')
        root.deiconify()
        return
    if not no_arg and result[0] is not None:
        messagebox.showinfo(message=f'Done. See {result[0].out} for details.')
    else:
        messagebox.showinfo(message='Done')
    try:
        window.withdraw()
    except tk.TclError:
        pass
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

        img_dir = resources.files(name) / 'data'
        photo_location1 = img_dir / 'button1.png'
        photo_location2 = img_dir / 'button2.png'
        photo_location3 = img_dir / 'button3.png'
        photo_location4 = img_dir / 'button4.png'
        global _img4
        _img4 = tk.PhotoImage(file=photo_location4)
        global _img0
        _img0 = tk.PhotoImage(file=photo_location1)
        global _img1
        _img1 = tk.PhotoImage(file=photo_location2)
        global _img2
        _img2 = tk.PhotoImage(file=photo_location3)

        top.geometry("800x450+400+0")
        move_to_center(top, 800, 450)
        top.minsize(120, 15)
        top.resizable(1, 1)
        top.title(name)
        top.configure(background="#edf0f3")
        top.configure(highlightbackground="#edf0f3")
        top.configure(highlightcolor="black")
        self.top = top

        self.help_b = my_button(self.top)
        self.help_b.place(relx=0.913, rely=0.067, height=40, width=40)
        self.help_b.configure(command=run_help)
        self.help_b.configure(image=_img4)

        self.gb2fasta_b = my_button(self.top)
        self.gb2fasta_b.place(relx=0.188, rely=0.288, height=100, width=100)
        self.gb2fasta_b.configure(command=ui_gb2fasta)
        self.gb2fasta_b.configure(image=_img0)

        self.evaluate_b = my_button(self.top)
        self.evaluate_b.place(relx=0.438, rely=0.288, height=100, width=100)
        self.evaluate_b.configure(command=ui_evaluate)
        self.evaluate_b.configure(image=_img1)

        self.primer_b = my_button(self.top)
        self.primer_b.place(relx=0.688, rely=0.288, height=100, width=100)
        self.primer_b.configure(command=ui_primer)
        self.primer_b.configure(image=_img2)

        self.visualize = my_button(self.top)
        self.visualize.place(relx=0.125, rely=0.865, height=30, width=250)
        self.visualize.configure(text='Visualize')
        self.visualize.configure(command=run_visualize(self, self.top))

        self.install_third_party = my_button(self.top)
        self.install_third_party.place(relx=0.613, rely=0.865, height=30,
                                       width=250)
        self.install_third_party.configure(text='Install third-party software')
        self.install_third_party.configure(command=run_install(self, self.top))

        self.gb2fasta_label = my_label(self.top, 14)
        self.gb2fasta_label.place(relx=0.188, rely=0.532, height=30, width=100)
        self.gb2fasta_label.configure(text='''GB2Fasta''')

        self.evaluate_label = my_label(self.top, 14)
        self.evaluate_label.place(relx=0.45, rely=0.532, height=30, width=100)
        self.evaluate_label.configure(text='''Evaluate''')

        self.primer_label = my_label(self.top, 14)
        self.primer_label.place(relx=0.713, rely=0.532, height=30, width=100)
        self.primer_label.configure(text='''Primer''')

        self.note_label = my_label(self.top, 14)
        self.note_label.place(relx=0.35, rely=0.643, height=35, width=300)
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
        self.style.map('.', background=[('selected', _compcolor),
                                        ('active', _ana2color)])

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

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.025, rely=0.013, relheight=0.46
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Input''')

        self.TSeparator1 = ttk.Separator(self.Labelframe1)
        self.TSeparator1.place(relx=0.524, rely=0.19, relheight=0.541
                               , bordermode='ignore')
        self.TSeparator1.configure(orient="vertical")

        self.gbfile_label = my_label(self.Labelframe1)
        self.gbfile_label.place(relx=0.050, rely=0.057, height=35, width=100,
                                bordermode='ignore')
        self.gbfile_label.configure(text='''Genbank files''')

        self.gb_entry = my_entry(self.Labelframe1)
        self.gb_entry.place(relx=0.241, rely=0.057, relheight=0.095,
                            relwidth=0.562, bordermode='ignore')
        self.gb_entry.configure(textvariable=self.gb)
        self.gb_entry_tooltip = ToolTip(self.gb_entry, 'gb format files')

        self.gb_file_b = my_button(self.Labelframe1, 12)
        self.gb_file_b.place(relx=0.82, rely=0.054, height=35, width=90
                             , bordermode='ignore')
        self.gb_file_b.configure(command=open_file(self.gb_entry, single=False))
        self.gb_file_b.configure(text='''Open''')

        self.gene_label = my_label(self.Labelframe1)
        self.gene_label.place(relx=0.035, rely=0.217, height=35, width=60,
                              bordermode='ignore')
        self.gene_label.configure(text='''Gene''')

        self.gene_entry = my_entry(self.Labelframe1, 'eg. rbcL')
        self.gene_entry.place(relx=0.201, rely=0.217, height=35,
                              relwidth=0.314, bordermode='ignore')
        self.gene_entry.configure(textvariable=self.gene)
        self.gene_entry_tooltip = ToolTip(self.gene_entry, 'gene name')

        self.taxon_label = my_label(self.Labelframe1)
        self.taxon_label.place(relx=0.035, rely=0.326, height=35, width=72,
                               bordermode='ignore')
        self.taxon_label.configure(text='Taxonomy')

        self.molecular_label = my_label(self.Labelframe1)
        self.molecular_label.place(relx=0.558, rely=0.217, height=35, width=80
                                   , bordermode='ignore')
        self.molecular_label.configure(text='''Molecular''')

        self.TCombobox_molecular = my_combobox(self.Labelframe1)
        self.TCombobox_molecular.place(relx=0.716, rely=0.217, height=35,
                                       width=150, bordermode='ignore')
        self.molecular_value_list = ['all', 'DNA', 'RNA', ]
        self.TCombobox_molecular.configure(values=self.molecular_value_list)
        self.TCombobox_molecular.configure(textvariable=self.molecular)
        self.TCombobox_molecular.current(0)

        self.group_label = my_label(self.Labelframe1)
        self.group_label.place(relx=0.558, rely=0.326, height=35, width=60
                               , bordermode='ignore')
        self.group_label.configure(text='''Group''')

        self.TCombobox_group = my_combobox(self.Labelframe1)
        self.TCombobox_group.place(relx=0.716, rely=0.326, height=35, width=150,
                                   bordermode='ignore')
        self.group_value_list = ['all', 'animals', 'plants', 'fungi',
                                 'protists',
                                 'bacteria', 'archaea', 'viruses', ]
        self.TCombobox_group.configure(values=self.group_value_list)
        self.TCombobox_group.configure(textvariable=self.group)
        self.TCombobox_group.current(0)

        self.organelle_label = my_label(self.Labelframe1)
        self.organelle_label.place(relx=0.035, rely=0.435, height=35, width=80
                                   , bordermode='ignore')
        self.organelle_label.configure(text='''Organelle''')

        self.TCombobox_og = my_combobox(self.Labelframe1)
        self.TCombobox_og.place(relx=0.201, rely=0.435, height=35, width=170,
                                bordermode='ignore')
        self.og_value_list = ['ignore', 'both', 'no', 'mitochondrion',
                              'chloroplast', ]
        self.TCombobox_og.configure(values=self.og_value_list)
        self.TCombobox_og.configure(textvariable=self.og)
        self.TCombobox_og.current(0)

        self.refseq_label = my_label(self.Labelframe1)
        self.refseq_label.place(relx=0.035, rely=0.541, height=35, width=70
                                , bordermode='ignore')
        self.refseq_label.configure(text='''RefSeq''')

        self.TCombobox_refseq = my_combobox(self.Labelframe1)
        self.TCombobox_refseq.place(relx=0.201, rely=0.541, height=35,
                                    width=150, bordermode='ignore')
        self.refseq_value_list = ['both', 'yes', 'no', ]
        self.TCombobox_refseq.configure(values=self.refseq_value_list)
        self.TCombobox_refseq.configure(textvariable=self.refseq)
        self.TCombobox_refseq.current(0)
        self.TCombobox_refseq_tooltip = ToolTip(self.TCombobox_refseq,
                                                'Use RefSeq or not')

        self.count_label = my_label(self.Labelframe1)
        self.count_label.place(relx=0.035, rely=0.649, height=35, width=60
                               , bordermode='ignore')
        self.count_label.configure(text='''Count''')

        self.count_entry = my_entry(self.Labelframe1, 'eg. 1000')
        self.count_entry.place(relx=0.201, rely=0.649, relheight=0.095
                               , relwidth=0.314, bordermode='ignore')
        self.count_entry.configure(textvariable=self.count)
        self.count_entry_tooltip = ToolTip(
            self.count_entry, 'numbers of sequences to download, 0 for no limit')

        self.len_label = my_label(self.Labelframe1)
        self.len_label.place(relx=0.558, rely=0.435, height=35, width=60
                             , bordermode='ignore')
        self.len_label.configure(text='''Length''')

        self.min_len_entry = my_entry(self.Labelframe1)
        self.min_len_entry.place(relx=0.686, rely=0.435, relheight=0.095
                                 , relwidth=0.087, bordermode='ignore')
        self.min_len_entry.configure(textvariable=self.min_len)
        self.min_len_entry_tooltip = ToolTip(self.min_len_entry,
                                             'sequence length limit')

        self.to_label = my_label(self.Labelframe1)
        self.to_label.place(relx=0.785, rely=0.435, height=35, width=36
                            , bordermode='ignore')
        self.to_label.configure(text='''to''')

        self.max_len_entry = my_entry(self.Labelframe1)
        self.max_len_entry.place(relx=0.839, rely=0.435, relheight=0.095
                                 , relwidth=0.14, bordermode='ignore')
        self.max_len_entry.configure(textvariable=self.max_len)
        self.min_len_entry.insert(0, '0')
        self.max_len_entry.insert(0, '300000')

        self.date_label = my_label(self.Labelframe1)
        self.date_label.place(relx=0.558, rely=0.541, height=35, width=60
                              , bordermode='ignore')
        self.date_label.configure(text='''Date''')

        self.date_start_entry = my_entry(self.Labelframe1)
        self.date_start_entry.place(relx=0.686, rely=0.541, relheight=0.095
                                    , relwidth=0.087, bordermode='ignore')
        self.date_start_entry.configure(textvariable=self.date_start)
        self.date_start_entry_tooltip = ToolTip(self.date_start_entry,
                                                'YYYY/MM/DD')

        self.to2_label = my_label(self.Labelframe1)
        self.to2_label.place(relx=0.789, rely=0.541, height=35, width=36
                             , bordermode='ignore')
        self.to2_label.configure(text='''to''')

        self.date_end_entry = my_entry(self.Labelframe1)
        self.date_end_entry.place(relx=0.839, rely=0.541, relheight=0.095
                                  , relwidth=0.136, bordermode='ignore')
        self.date_end_entry.configure(textvariable=self.date_end)
        self.date_end_entry_tooltip = ToolTip(self.date_end_entry, 'YYYY/MM/DD')

        self.exclude_label = my_label(self.Labelframe1)
        self.exclude_label.place(relx=0.558, rely=0.649, height=35, width=60
                                 , bordermode='ignore')
        self.exclude_label.configure(text='''Exclude''')

        self.exclude_entry = my_entry(self.Labelframe1, 'eg. Fungi')
        self.exclude_entry.place(relx=0.686, rely=0.649, relheight=0.095
                                 , relwidth=0.279, bordermode='ignore')
        self.exclude_entry.configure(textvariable=self.exclude)
        self.exclude_entry_tooltip = ToolTip(self.exclude_entry,
                                             'excluded taxon, can be empty')

        self.query_label = my_label(self.Labelframe1)
        self.query_label.place(relx=0.035, rely=0.813, height=35, width=60
                               , bordermode='ignore')
        self.query_label.configure(text='''Query''')

        self.query_entry = my_entry(self.Labelframe1,
                                    'eg. rpoB[gene] AND txid9721[organism]')
        self.query_entry.place(relx=0.201, rely=0.813, relheight=0.095
                               , relwidth=0.785, bordermode='ignore')
        self.query_entry.configure(textvariable=self.query)
        self.query_entry_tooltip = ToolTip(self.query_entry, 'Entrez query string')

        self.taxon_entry = my_entry(self.Labelframe1, 'eg. Fabaceae')
        self.taxon_entry.place(relx=0.201, rely=0.326, relheight=0.095
                               , relwidth=0.316, bordermode='ignore')
        self.taxon_entry.configure(textvariable=self.taxon)
        self.taxon_entry_tooltip = ToolTip(
            self.taxon_entry,
            'Taxonomy name (Zea mays, Fabaceae, etc.) or Txid (eg. txid9721)')

        self.out_label = my_label(self.top)
        self.out_label.place(relx=0.05, rely=0.563, height=36
                             , width=60)
        self.out_label.configure(text='''Output''')

        self.out_entry = my_entry(self.top, 'eg. F:/mywork2/result')
        self.out_entry.place(relx=0.167, rely=0.563, relheight=0.045
                             , relwidth=0.633)
        self.out_entry.configure(textvariable=self.out)
        self.out_entry_tooltip = ToolTip(self.out_entry, 'output folder')
        self.out_b = my_button(self.top, 12)
        self.out_b.place(relx=0.833, rely=0.563, height=35, width=80)
        self.out_b.configure(command=open_file(self.out_entry, type_='folder'))
        self.out_b.configure(text='''Open''')

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.025, rely=0.65, relheight=0.176
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Advance''')

        self.expand_label = my_label(self.Labelframe1)
        self.expand_label.place(relx=0.185, rely=0.142, height=35
                                , width=60, bordermode='ignore')
        self.expand_label.configure(text='''Expand''')

        self.expand_entry = my_entry(self.Labelframe1, '0')
        self.expand_entry.place(relx=0.297, rely=0.142, relheight=0.248
                                , relwidth=0.187, bordermode='ignore')
        self.expand_entry.configure(textvariable=self.expand)
        self.expand_entry_tooltip = ToolTip(
            self.expand_entry,
            'expand how many bp to upstream/downstream for primer design')

        self.max_name_len_label = my_label(self.Labelframe1)
        self.max_name_len_label.place(relx=0.08, rely=0.426, height=35
                                      , width=130, bordermode='ignore')
        self.max_name_len_label.configure(text='''Max name length''')

        self.max_frag_len_label = my_label(self.Labelframe1)
        self.max_frag_len_label.place(relx=0.045, rely=0.709, height=35
                                      , width=150, bordermode='ignore')
        self.max_frag_len_label.configure(text='''Max fragment length''')

        self.max_name_len_entry = my_entry(self.Labelframe1)
        self.max_name_len_entry.place(relx=0.297, rely=0.426, relheight=0.248
                                      , relwidth=0.187, bordermode='ignore')
        self.max_name_len_entry.configure(textvariable=self.max_name_len)
        self.max_name_len_entry_tooltip = ToolTip(self.max_name_len_entry,
                                                  'max feature name length')

        self.max_gene_len_entry = my_entry(self.Labelframe1)
        self.max_gene_len_entry.place(relx=0.297, rely=0.709, relheight=0.248
                                      , relwidth=0.187, bordermode='ignore')
        self.max_gene_len_entry.configure(textvariable=self.max_gene_len)
        self.max_gene_len_entry_tooltip = ToolTip(self.max_gene_len_entry,
                                                  'max fragment sequence length')
        self.max_name_len_entry.insert(0, '100')
        self.max_gene_len_entry.insert(0, '20000')

        self.TSeparator2 = ttk.Separator(self.Labelframe1)
        self.TSeparator2.place(relx=0.524, rely=0.135, relheight=0.78
                               , bordermode='ignore')
        self.TSeparator2.configure(orient="vertical")

        self.allow_repeat_button = my_checkbutton(self.Labelframe1)
        self.allow_repeat_button.place(relx=0.576, rely=0.142, relheight=0.248
                                       , relwidth=0.276, bordermode='ignore')
        self.allow_repeat_button.configure(text='''Allow repeat''')
        self.allow_repeat_button.configure(variable=self.allow_repeat)

        self.allow_invert_repeat_button = my_checkbutton(self.Labelframe1)
        self.allow_invert_repeat_button.place(relx=0.576, rely=0.709,
                                              relheight=0.248
                                              , relwidth=0.361,
                                              bordermode='ignore')
        self.allow_invert_repeat_button.configure(
            text='''Allow invert repeat''')
        self.allow_invert_repeat_button.configure(
            variable=self.allow_invert_repeat)

        self.allow_mosaic_button = my_checkbutton(self.Labelframe1)
        self.allow_mosaic_button.place(relx=0.576, rely=0.426, relheight=0.248
                                       , relwidth=0.361, bordermode='ignore')
        self.allow_mosaic_button.configure(text='''Allow mosaic repeat''')
        self.allow_mosaic_button.configure(variable=self.allow_mosaic_repeat)

        self.no_divide_b = my_checkbutton(self.top)
        self.no_divide_b.place(relx=0.733, rely=0.488, relheight=0.044
                               , relwidth=0.215)
        self.no_divide_b.configure(text='''No divide''')
        self.no_divide_b.configure(variable=self.no_divide)

        self.rename_b = my_checkbutton(self.top)
        self.rename_b.place(relx=0.467, rely=0.488, relheight=0.044
                            , relwidth=0.243)
        self.rename_b.configure(text='''Rename gene''')
        self.rename_b.configure(variable=self.rename)

        self.unique_label = my_label(self.top)
        self.unique_label.place(relx=0.05, rely=0.488, height=35
                                , width=69)
        self.unique_label.configure(text='''Unique''')

        self.TCombobox_unique = my_combobox(self.top)
        self.TCombobox_unique.place(relx=0.167, rely=0.488, relheight=0.044
                                    , relwidth=0.245)
        self.unique_value_list = ['longest', 'first', 'no', ]
        self.TCombobox_unique.configure(values=self.unique_value_list)
        self.TCombobox_unique.configure(textvariable=self.unique)
        self.TCombobox_unique.current(0)
        self.TCombobox_unique_tooltip = ToolTip(
            self.TCombobox_unique, 'methods to remove redundant records')

        self.run_b = my_button(self.top, 14)
        self.run_b.place(relx=0.333, rely=0.863, height=40, width=189)
        self.run_b.configure(command=run_gb2fasta(self, self.top))
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

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.0, relheight=0.333
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Input''')

        self.unalign_label = my_label(self.Labelframe1)
        self.unalign_label.place(relx=0.03, rely=0.133, height=35, width=160
                                 , bordermode='ignore')
        self.unalign_label.configure(text='''Unaligned FASTA files''')
        self.TLabel1_tooltip = ToolTip(self.unalign_label, 'unaligned')

        self.fasta_entry = my_entry(self.Labelframe1, 'eg. c:/mydata/a.fasta')
        self.fasta_entry.place(relx=0.314, rely=0.133, relheight=0.233
                               , relwidth=0.489, bordermode='ignore')
        self.fasta_entry.configure(textvariable=self.fasta)
        self.fasta_entry_tooltip = ToolTip(self.fasta_entry,
                                           'one or more unaligned fasta files')

        self.open_btn = my_button(self.Labelframe1)
        self.open_btn.place(relx=0.82, rely=0.133, height=35, width=90
                            , bordermode='ignore')
        self.open_btn.configure(command=open_file(self.fasta_entry,
                                                  single=False))
        self.open_btn.configure(text='''Open''')

        self.unalign_label2 = my_label(self.Labelframe1)
        self.unalign_label2.place(relx=0.03, rely=0.4, height=35, width=160
                                  , bordermode='ignore')
        self.unalign_label2.configure(text='''Unaligned FASTA folder''')
        self.TLabel1_tooltip = ToolTip(self.unalign_label2, 'unaligned')

        self.fasta_folder_entry = my_entry(self.Labelframe1,
                                           'eg. d:/Result/Divide')
        self.fasta_folder_entry.place(relx=0.314, rely=0.4, relheight=0.233
                                      , relwidth=0.489, bordermode='ignore')
        self.fasta_folder_entry.configure(textvariable=self.fasta_folder)
        self.fasta_folder_entry_tooltip = ToolTip(
            self.fasta_folder_entry,
            'a folder containing one or more unaligned fasta files')
        self.open1_btn = my_button(self.Labelframe1)
        self.open1_btn.place(relx=0.82, rely=0.4, height=35, width=90
                             , bordermode='ignore')
        self.open1_btn.configure(command=open_file(self.fasta_folder_entry,
                                                   type_='folder'))
        self.open1_btn.configure(text='''Open''')

        self.align_label = my_label(self.Labelframe1)
        self.align_label.place(relx=0.03, rely=0.667, height=35, width=150
                               , bordermode='ignore')
        self.align_label.configure(text='''Aligned FASTA files''')
        self.TLabel1_3_tooltip = ToolTip(self.align_label, 'aligned')

        self.aln_entry = my_entry(self.Labelframe1, 'eg. e:/Result/Aligned')
        self.aln_entry.place(relx=0.314, rely=0.667, relheight=0.233
                             , relwidth=0.489, bordermode='ignore')
        self.aln_entry.configure(textvariable=self.aln)
        self.aln_entry_tooltip = ToolTip(self.aln_entry,
                                         'one or more aligned fasta files')
        self.open2_btn = my_button(self.Labelframe1)
        self.open2_btn.place(relx=0.82, rely=0.667, height=35, width=90
                             , bordermode='ignore')
        self.open2_btn.configure(command=open_file(self.aln_entry,
                                                   single=False))
        self.open2_btn.configure(text='''Open''')

        self.out_label = my_label(self.top)
        self.out_label.place(relx=0.117, rely=0.356, height=35, width=100)
        self.out_label.configure(text='''Output folder''')
        self.TLabel1_3_1_tooltip = ToolTip(self.out_label, 'output')

        self.out_entry = my_entry(self.top, 'eg. H:/project9/Result')
        self.out_entry.place(relx=0.317, rely=0.356, relheight=0.078
                             , relwidth=0.467)
        self.out_entry.configure(textvariable=self.out)
        self.out_entry_tooltip = ToolTip(self.out_entry, 'unaligned fasta files')

        self.open3_btn = my_button(self.top)
        self.open3_btn.place(relx=0.8, rely=0.356, height=35, width=90)
        self.open3_btn.configure(command=open_file(self.out_entry,
                                                   type_='folder'))
        self.open3_btn.configure(text='''Open''')

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.467, relheight=0.156
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Sliding window''')

        self.Checkbutton1 = my_checkbutton(self.Labelframe1)
        self.Checkbutton1.place(relx=0.005, rely=0.429, relheight=0.3,
                                relwidth=0.262, bordermode='ignore')
        self.Checkbutton1.configure(text='''Skip sliding window''')
        self.Checkbutton1.configure(variable=self.quick)

        self.window_size_label = my_label(self.Labelframe1)
        self.window_size_label.place(relx=0.290, rely=0.429, height=22,
                                     width=90, bordermode='ignore')
        self.window_size_label.configure(text='Window size')

        self.step_len_label = my_label(self.Labelframe1)
        self.step_len_label.place(relx=0.681, rely=0.429, height=22, width=80
                                  , bordermode='ignore')
        self.step_len_label.configure(text='''Step length''')

        self.size_entry = my_entry(self.Labelframe1)
        self.size_entry.place(relx=0.489, rely=0.357, relheight=0.5
                              , relwidth=0.14, bordermode='ignore')
        self.size_entry.configure(textvariable=self.size)

        self.step_entry = my_entry(self.Labelframe1)
        self.step_entry.place(relx=0.838, rely=0.357, relheight=0.5
                              , relwidth=0.14, bordermode='ignore')
        self.step_entry.configure(textvariable=self.step)
        self.size_entry.insert(0, '500')
        self.step_entry.insert(0, '50')

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.644, relheight=0.156
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Advance''')

        self.Checkbutton1 = my_checkbutton(self.Labelframe1)
        self.Checkbutton1.place(relx=0.087, rely=0.429, relheight=0.3
                                , relwidth=0.354, bordermode='ignore')
        self.Checkbutton1.configure(text='''Ignore gaps in alignment''')
        self.Checkbutton1.configure(variable=self.ig)

        self.Checkbutton1_2 = my_checkbutton(self.Labelframe1)
        self.Checkbutton1_2.place(relx=0.593, rely=0.429, relheight=0.3
                                  , relwidth=0.354, bordermode='ignore')
        self.Checkbutton1_2.configure(text='''Ignore ambiguous bases''')
        self.Checkbutton1_2.configure(variable=self.iab)

        self.run_b = my_button(self.top, 14)
        self.run_b.place(relx=0.333, rely=0.867, height=40, width=189)
        self.run_b.configure(command=run_evaluate(self, self.top))
        self.run_b.configure(text='''Run''')


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
        self.style.map('.', background=[('selected', _compcolor),
                                        ('active', _ana2color)])

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

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.02, relheight=0.2
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Input''')

        self.aln_fasta_label = my_label(self.Labelframe1)
        self.aln_fasta_label.place(relx=0.052, rely=0.2, height=35, width=130
                                   , bordermode='ignore')
        self.aln_fasta_label.configure(text='''Aligned FASTA files''')
        self.TLabel1_tooltip = ToolTip(self.aln_fasta_label, 'unaligned')

        self.aln_entry = my_entry(self.Labelframe1, 'eg. f:/aligned/atpB.aln')
        self.aln_entry.place(relx=0.314, rely=0.2, relheight=0.35,
                             relwidth=0.489
                             , bordermode='ignore')
        self.aln_entry.configure(textvariable=self.aln)
        self.aln_entry_tooltip = ToolTip(self.aln_entry,
                                         'one or more unaligned fasta files')

        self.out_b = my_button(self.Labelframe1)
        self.out_b.place(relx=0.82, rely=0.2, height=35, width=90
                         , bordermode='ignore')
        self.out_b.configure(command=open_file(self.aln_entry, single=False))
        self.out_b.configure(text='''Open''')

        self.aln_folder_label = my_label(self.Labelframe1)
        self.aln_folder_label.place(relx=0.052, rely=0.6, height=35, width=160
                                    , bordermode='ignore')
        self.aln_folder_label.configure(text='''Aligned FASTA folder''')
        self.TLabel1_tooltip = ToolTip(self.aln_folder_label,
                                       'Folder with aligned fasta files')

        self.aln_folder_entry = my_entry(self.Labelframe1,
                                         'eg. g:/mywork3/AlignedResult')
        self.aln_folder_entry.place(relx=0.314, rely=0.6, relheight=0.35
                                    , relwidth=0.489, bordermode='ignore')
        self.aln_folder_entry.configure(textvariable=self.aln_folder)
        self.aln_folder_entry_tooltip = ToolTip(
            self.aln_folder_entry,
            'a folder containing one or more unaligned fasta files')

        self.folder_b = my_button(self.Labelframe1)
        self.folder_b.place(relx=0.82, rely=0.6, height=35, width=90
                            , bordermode='ignore')
        self.folder_b.configure(command=open_file(self.aln_folder_entry,
                                                  type_='folder'))
        self.folder_b.configure(text='''Open''')

        self.out_label = my_label(self.top)
        self.out_label.place(relx=0.07, rely=0.24, height=35, width=150)
        self.out_label.configure(text='''Output folder''')
        self.TLabel1_3_1_tooltip = ToolTip(self.out_label, 'Output')

        self.out_entry = my_entry(self.top)
        self.out_entry.place(relx=0.317, rely=0.24, relheight=0.07
                             , relwidth=0.467)
        self.out_entry.configure(textvariable=self.out)

        self.out_b = my_button(self.top)
        self.out_b.place(relx=0.8, rely=0.24, height=35, width=90)
        self.out_b.configure(command=open_file(self.out_entry, type_='folder'))
        self.out_b.configure(text='''Open''')

        self.Labelframe1 = my_labelframe(self.top)
        self.Labelframe1.place(relx=0.017, rely=0.34, relheight=0.38
                               , relwidth=0.955)
        self.Labelframe1.configure(text='''Advance''')

        self.coverage_label = my_label(self.Labelframe1)
        self.coverage_label.place(relx=0.07, rely=0.158, height=35
                                  , width=70, bordermode='ignore')
        self.coverage_label.configure(text='''Coverage''')

        self.coverage_entry = my_entry(self.Labelframe1)
        self.coverage_entry.place(relx=0.297, rely=0.142, relheight=0.184
                                  , relwidth=0.192, bordermode='ignore')
        self.coverage_entry.configure(textvariable=self.coverage)
        self.coverage_entry_tooltip = ToolTip(
            self.coverage_entry, 'minimal coverage of primer on alignment')
        self.coverage_entry.insert(0, '0.5')

        self.mismatch_label = my_label(self.Labelframe1)
        self.mismatch_label.place(relx=0.07, rely=0.368, height=35
                                  , width=120, bordermode='ignore')
        self.mismatch_label.configure(text='''Mismatch''')

        self.res_label = my_label(self.Labelframe1)
        self.res_label.place(relx=0.07, rely=0.579, height=35
                             , width=130, bordermode='ignore')
        self.res_label.configure(text='''Resolution''')

        self.mismatch_entry = my_entry(self.Labelframe1)
        self.mismatch_entry.place(relx=0.297, rely=0.368, relheight=0.184
                                  , relwidth=0.192, bordermode='ignore')
        self.mismatch_entry.configure(textvariable=self.mismatch)
        self.mismatch_entry_tooltip = ToolTip(
            self.mismatch_entry, 'maximum mismatch bases in primer')

        self.resolution_entry = my_entry(self.Labelframe1)
        self.resolution_entry.place(relx=0.297, rely=0.579, relheight=0.184
                                    , relwidth=0.192, bordermode='ignore')
        self.resolution_entry.configure(textvariable=self.resolution)
        self.resolution_entry_tooltip = ToolTip(
            self.resolution_entry, 'minimal resolution of amplified fragment')
        self.mismatch_entry.insert(0, '4')
        self.resolution_entry.insert(0, '0.3')

        self.TSeparator2 = ttk.Separator(self.Labelframe1)
        self.TSeparator2.place(relx=0.524, rely=0.137, relheight=0.779
                               , bordermode='ignore')
        self.TSeparator2.configure(orient="vertical")

        self.topn_label = my_label(self.Labelframe1)
        self.topn_label.place(relx=0.07, rely=0.789, height=35
                              , width=130, bordermode='ignore')
        self.topn_label.configure(text='''Top n''')

        self.top_n_entry = my_entry(self.Labelframe1)
        self.top_n_entry.place(relx=0.297, rely=0.789, relheight=0.184
                               , relwidth=0.192, bordermode='ignore')
        self.top_n_entry.configure(textvariable=self.top_n)
        self.top_n_entry_tooltip = ToolTip(self.top_n_entry,
                                           'Only keep top best primers')
        self.top_n_entry.insert(0, '1')

        self.primer_len_label = my_label(self.Labelframe1)
        self.primer_len_label.place(relx=0.541, rely=0.158, height=35
                                    , width=100, bordermode='ignore')
        self.primer_len_label.configure(text='''Primer length''')

        self.pmin_entry = my_entry(self.Labelframe1)
        self.pmin_entry.place(relx=0.716, rely=0.158, relheight=0.184
                              , relwidth=0.087, bordermode='ignore')
        self.pmin_entry.configure(textvariable=self.pmin)

        self.to_label1 = my_label(self.Labelframe1)
        self.to_label1.place(relx=0.803, rely=0.158, height=35, width=36
                             , bordermode='ignore')
        self.to_label1.configure(text='''to''')

        self.pmax_entry = my_entry(self.Labelframe1)
        self.pmax_entry.place(relx=0.855, rely=0.158, relheight=0.184
                              , relwidth=0.122, bordermode='ignore')
        self.pmax_entry.configure(textvariable=self.pmax)
        self.pmin_entry.insert(0, '20')
        self.pmax_entry.insert(0, '30')

        self.amp_len_label = my_label(self.Labelframe1)
        self.amp_len_label.place(relx=0.541, rely=0.368, height=35
                                 , width=120, bordermode='ignore')
        self.amp_len_label.configure(text='''Amplicon size''')

        self.amin_entry = my_entry(self.Labelframe1)
        self.amin_entry.place(relx=0.716, rely=0.368, relheight=0.184
                              , relwidth=0.087, bordermode='ignore')
        self.amin_entry.configure(textvariable=self.amin)
        self.amin_entry_tooltip = ToolTip(self.amin_entry,
                                          'including primer length')

        self.to2_label = my_label(self.Labelframe1)
        self.to2_label.place(relx=0.803, rely=0.368, height=35, width=36
                             , bordermode='ignore')
        self.to2_label.configure(text='''to''')

        self.amax_entry = my_entry(self.Labelframe1)
        self.amax_entry.place(relx=0.855, rely=0.368, relheight=0.184
                              , relwidth=0.122, bordermode='ignore')
        self.amax_entry.configure(textvariable=self.amax)
        self.amin_entry.insert(0, '300')
        self.amax_entry.insert(0, '800')

        self.sliding_window_label = my_label(self.Labelframe1)
        self.sliding_window_label.place(relx=0.541, rely=0.579, height=35,
                                        width=130
                                        , bordermode='ignore')
        self.sliding_window_label.configure(text='''Sliding window size''')

        self.sliding_window2_label = my_label(self.Labelframe1)
        self.sliding_window2_label.place(relx=0.541, rely=0.789, height=35,
                                         width=130
                                         , bordermode='ignore')
        self.sliding_window2_label.configure(text='''Sliding window step''')

        self.size_entry = my_entry(self.Labelframe1)
        self.size_entry.place(relx=0.803, rely=0.579, relheight=0.184
                              , relwidth=0.14, bordermode='ignore')
        self.size_entry.configure(textvariable=self.size)

        self.step_entry = my_entry(self.Labelframe1)
        self.step_entry.place(relx=0.803, rely=0.789, relheight=0.184
                              , relwidth=0.14, bordermode='ignore')
        self.step_entry.configure(textvariable=self.step)
        self.size_entry.insert(0, '500')
        self.step_entry.insert(0, '50')

        self.run_b = my_button(self.top, 14)
        self.run_b.place(relx=0.367, rely=0.82, height=40, width=180)
        self.run_b.configure(command=run_primer(self, self.top))
        self.run_b.configure(text='''Run''')


class ToolTip(tk.Toplevel):
    """ Provides a ToolTip widget for Tkinter. """

    def __init__(self, wdgt, msg=None, tooltip_font='TkDefaultFont 12',
                 msgFunc=None, delay=0.01, follow=True):
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
    value = value.get()
    # remove placeholder
    if isinstance(value, str) and value.startswith('eg. '):
        value = ''
    if value:
        if is_bool:
            value_str = f'{name} '
        else:
            value_str = f'{name} {value} '
    arg_str += value_str
    log.debug(f'{name} {value_str}')
    return arg_str


def run_gb2fasta(win, t: tk.Toplevel):
    # todo: test options and functions
    def f():
        nonlocal win
        arg_str = ''
        arg_str = get_arg_str(win.gb, '-gb', arg_str)
        arg_str = get_arg_str(win.gene, '-gene', arg_str)
        arg_str = get_arg_str(win.molecular, '-molecular', arg_str)
        arg_str = get_arg_str(win.group, '-group', arg_str)
        arg_str = get_arg_str(win.og, '-og', arg_str)
        arg_str = get_arg_str(win.refseq, '-refseq', arg_str)
        arg_str = get_arg_str(win.count, '-count', arg_str)
        arg_str = get_arg_str(win.min_len, '-min_len', arg_str)
        arg_str = get_arg_str(win.max_len, '-max_len', arg_str)
        arg_str = get_arg_str(win.date_start, '-date_start', arg_str)
        arg_str = get_arg_str(win.date_end, '-date_end', arg_str)
        arg_str = get_arg_str(win.exclude, '-exclude', arg_str)
        arg_str = get_arg_str(win.query, '-query', arg_str)
        arg_str = get_arg_str(win.taxon, '-taxon', arg_str)
        arg_str = get_arg_str(win.out, '-out', arg_str)
        arg_str = get_arg_str(win.expand, '-expand', arg_str)
        arg_str = get_arg_str(win.max_name_len, '-max_name_len', arg_str)
        arg_str = get_arg_str(win.max_gene_len, '-max_gene_len', arg_str)
        arg_str = get_arg_str(win.unique, '-unique', arg_str)
        arg_str = get_arg_str(win.allow_repeat, '-allow_repeat', arg_str,
                              is_bool=True)
        arg_str = get_arg_str(win.allow_invert_repeat, '-allow_invert_repeat',
                              arg_str,
                              is_bool=True)
        arg_str = get_arg_str(win.allow_mosaic_repeat, '-allow_mosaic_repeat',
                              arg_str,
                              is_bool=True)
        arg_str = get_arg_str(win.no_divide, '-no_divide', arg_str,
                              is_bool=True)
        arg_str = get_arg_str(win.rename, '-rename', arg_str, is_bool=True)
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s * 0.618)}+{w // 3}+{h // 3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(gb2fasta_main, arg_str, run),
                             daemon=True)
        r.start()

    return f


def run_evaluate(win, t: tk.Toplevel):
    # todo: test options and functions
    def f():
        nonlocal win
        arg_str = ''
        arg_str = get_arg_str(win.fasta, '-fasta', arg_str)
        arg_str = get_arg_str(win.fasta_folder, '-fasta_folder', arg_str)
        arg_str = get_arg_str(win.aln, '-aln', arg_str)
        arg_str = get_arg_str(win.out, '-out', arg_str)
        arg_str = get_arg_str(win.size, '-size', arg_str)
        arg_str = get_arg_str(win.step, '-step', arg_str)
        arg_str = get_arg_str(win.quick, '-quick', arg_str, is_bool=True)
        arg_str = get_arg_str(win.ig, '-ig', arg_str, is_bool=True)
        arg_str = get_arg_str(win.iab, '-iab', arg_str, is_bool=True)
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s * 0.618)}+{w // 3}+{h // 3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(evaluate_main, arg_str, run),
                             daemon=True)
        r.start()

    return f


def run_primer(win, t: tk.Toplevel):
    # todo: test options and functions
    def f():
        nonlocal win
        arg_str = ''
        arg_str = get_arg_str(win.aln, '-aln', arg_str)
        arg_str = get_arg_str(win.aln_folder, '-aln_folder', arg_str)
        arg_str = get_arg_str(win.out, '-out', arg_str)
        arg_str = get_arg_str(win.coverage, '-coverage', arg_str)
        arg_str = get_arg_str(win.mismatch, '-mismatch', arg_str)
        arg_str = get_arg_str(win.resolution, '-resolution', arg_str)
        arg_str = get_arg_str(win.top_n, '-top_n', arg_str)
        arg_str = get_arg_str(win.pmin, '-pmin', arg_str)
        arg_str = get_arg_str(win.pmax, '-pmax', arg_str)
        arg_str = get_arg_str(win.amin, '-amin', arg_str)
        arg_str = get_arg_str(win.amax, '-amax', arg_str)
        arg_str = get_arg_str(win.size, '-size', arg_str)
        arg_str = get_arg_str(win.step, '-step', arg_str)
        arg_str += '-primer'
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s * 0.618)}+{w // 3}+{h // 3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(primer_main, arg_str, run),
                             daemon=True)
        r.start()

    return f


def run_help():
    url = 'https://github.com/wpwupingwp/OGU'
    webbrowser.open(url, new=2)


def func_wrap(t: tk.Toplevel, func):
    def _f():
        nonlocal t, func
        t.withdraw()
        w, h = root.winfo_screenwidth(), root.winfo_screenheight()
        s = min(w, h) // 2
        size = f'{s}x{int(s * 0.618)}+{w // 3}+{h // 3}'
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap, args=(func, '', run, True),
                             daemon=True)
        r.start()
    return _f


def run_visualize(win, t: tk.Toplevel):
    def v():
        # with (resources.files(name) / 'data') as v_folder:
        v_folder = resources.files(name) / 'data'
        system = platform.system()
        if system == "Windows":
            os.startfile(v_folder)
        elif system == "Darwin":  # macOS
            os.system(f"open '{v_folder}'")
        elif system == "Linux":
            os.system(f"xdg-open '{v_folder}'")
        else:
            log.error('Unsupported system')
    v_func = func_wrap(t, v)
    return v_func


def run_install(win, t: tk.Toplevel):
    i_func = func_wrap(t, get_all_third_party)
    return i_func


def ui_main():
    global root
    root = tk.Tk()
    root.protocol('WM_DELETE_WINDOW', root.destroy)
    root.title('Organelle Genome Utilities')
    # Creates a toplevel widget.
    global _top1, _w1
    _top1 = root
    _w1 = Root(_top1)
    root.mainloop()


if __name__ == '__main__':
    ui_main()