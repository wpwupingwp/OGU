#!/usr/bin/python3

with open('/usr/share/dict/words', 'r') as In:
    wordlist = In.read().split(sep='\n')
vowel = set('aeiouAEIOU')
all_vowel = list()

def is_all_vowel(word):
    for letter in word:
        if letter not in vowel:
            return False
    return True

for word in wordlist:
    if is_all_vowel(word) is True:
        all_vowel.append(word)

sorted(all_vowel, key = lambda x:len(x))

print('word length')
for word in all_vowel:
    print(word,len(word))
