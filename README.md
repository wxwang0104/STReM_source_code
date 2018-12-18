# Summary
This repository contains the STReM source code for a simple demonstration of 3D decomposition

# Copyright
Copyright (c) 2018, Wenxiao Wang All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in
     the documentation and/or other materials provided with the distribution

# STReM_source_code
Software:   Matlab 2009 or higher version

Run 
```
demo_cos.m
```
for a simple demo

# Simulate a cosin shaped trajectory 

1. The images in the left show the PSFs of emitters at different positions; 
      The images in the right show the cumulated images
![Raw img](processes.png)

2. The final captured image that is compressed with time information
![Raw img](raw.png)  
  
# Deconvolve the raw trajectory with 3D DH PSFs

The results are overlaid the recovered trajectory with raw image. 

![Raw img](recovery.png)
