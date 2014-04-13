#Guess the number, need simplegui
import simplegui
import random

secret=0
range=0
guesses=0
guessed=0

def new_game(range):
    global secret
    secret=random.randrange(0,range)
    global guessed
    guessed=0
    print "My secret number is an integer that smaller than %d. Guess it!\n" %range

def range100():
    global guesses
    guesses=7
    global range
    range=100
    new_game(range)

def range1000():
    global guesses
    guesses=10
    global range
    range=1000
    new_game(range)
    
def get_text(text_input):
    global guessed
    guess=int(text_input)
    if guess>=0 and guess<range:
        guessed=guessed+1
        input_guess(guess)
    else:
        guessed=guessed+1
        print "Seriously? You just waste one guess! Give a correct input and you have only %d chances now." %(guesses-guessed)

def input_guess(guess):
    if guesses==guessed:
        print "Guesses used up. You lose.\n"
        new_game(range)
        return 
    elif guessed==guesses-1:
        print "Last chance, be careful!"
    else:
        print "Now you have %d guesses." %(guesses-guessed)

    if guess==secret:
        print "%d! After %d times tries, you win!\n" %(guess,guessed)
        new_game(range)
    elif guess<secret:
        print "%d? You need more courage! Guess a bigger number." %guess
    else:
        print "%d? You are reckless this time. Please try a smaller one." %guess
    
frame=simplegui.create_frame("Guess the number",200,200)

button1=frame.add_button("Range [0,100)",range100,50)
button2=frame.add_button("Range [0,1000)",range1000,70)
inp=frame.add_input("Enter the range:",get_text,50)

range100()
frame.start
