#!/usr/bin/python3

with open('/usr/share/dict/words', 'r') as In:
    wordlist = In.read().split(sep='\n')
vowel = set('aeiouAEIOU')
under = set('underUNDER')
all_vowel = list()
all_under = list()


def is_all_vowel(word, letters):
    for letter in word:
        if letter not in letters:
            return False
    return True


for word in wordlist:
    if is_all_vowel(word, vowel) is True:
        all_vowel.append(word)
    if is_all_vowel(word, under) is True:
        all_under.append(word)

all_under.sort(key=lambda x: len(x))
print('word length')
for word in all_under:
    print(word, len(word))
