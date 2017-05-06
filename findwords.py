#!/usr/bin/python3

with open('/usr/share/dict/words', 'r') as In:
    wordlist = In.read().split(sep='\n')
vowel = set('aeiouAEIOU')
all_vowel = list()

def is_all_vowel(word, letters):
    for letter in word:
        if letter not in letters:
            return False
    return True


for word in wordlist:
    if is_all_vowel(word, vowel) is True:
        all_vowel.append(word)

all_vowel.sort(key=lambda x: len(x))
print('word\tlength')
for word in all_vowel:
    print(word, len(word))
