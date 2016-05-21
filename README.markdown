#Differential ChipSeq Peak Finder

This program helps users analyze differential expression from ChipSeq data. Different Peaks from both the files are located and annotated with relevant gene, promoter, and enhancer info. Currently only mouse ChipSeq files are supported. The input sequence files need to be sorted and output files will be in a `.xls` format.

##Demo

[![Differential ChipSeq Peak Finder Video](http://img.youtube.com/vi/5XIOwX2ETU8/0.jpg)](https://www.youtube.com/watch?v=5XIOwX2ETU8&feature=youtu.be)

Note: The newer program allows you to adjust the tolerance (maximum distance between the centers of 2 peaks for it to still be considered 1 peak)

##Installation 

The program needs to have python version 2 installed. This program will not work with python version 3.

You need to first install the required python packages
```
sudo pip install -r requirements.txt
```

##Running

You can view the help dialog anytime by typing in:

```
python ChipSeqAlg.py -help 
```

There are 3 different way you can run this program 

###A short(1 session) web server

This server will be terminated as soon as you exit the terminal

Start the server by typing:

```
python ChipSeqAlg.py --website
``` 

The program will then output the url it is running on in the terminal. The program will run on localhost, so it can only be accessed by your computer.

###Using the Terminal

You can run the analysis manually by typing:

```
python ChipSeqAlg.py [ChipSeqFile1] [ChipSeqFile2]
```

###Using a long term web server

If you own a web server, you can run this program on it. First start supervisord: 

```
sudo supervisord -c supervisord.conf 
```
Note: If you use `python2` instead of `python` you edit this in the supervisord.conf by changing line 24.

The program should be running on port `5000`
Now visit the program at `[url of server]:5000`

You can check the status of the program by running:

```
sudo supervisorctl status ChipSeqAlg
```

You can stop the program by running:

```
sudo supervisorctl stop ChipSeqAlg
```

You can start the program by running:
```
sudo supervisorctl start ChipSeqAlg
```

You can restart the program by running:
```
sudo supervisorctl restart ChipSeqAlg
```