#!/usr/bin/python
import sys, os
import json 
from awesome_print import ap 
handles = sys.argv[1:]
print handles

json_rules = {}
output = {}
lines = {}

for handle in handles: 
  file = open(handle)
  j = json.load(file)
  rules = j['rules'] 
  for rule in rules: 
    output[rule] = {}
    lines[rule]=[]

for handle in handles:
  file = open(handle)
  j = json.load(file)
  rules = j['rules']
  for rule in rules:
    output[rule][handle]=j['rules'][rule]

for rule in output: 
  for handle in output[rule]:
    time = output[rule][handle]['mean-runtime'] 
    lines[rule].append(repr(time))

outlines = [['rule']]
for handle in handles:
  outlines[0].append(handle)

for line in lines:
  outline = [line] + lines[line] 
  outlines.append(outline)

for line in outlines: 
  print('\t'.join(line))
