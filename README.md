# BME 160 Code Repository
This is a repository of the code labs I have done in BME 160 - Research Programming in the Life Sciences at UCSC. Below is an excerpt of the syllabus describing the lab assignments done in this class,
the learning objectives for this class, and information on how the format of the class with a brief description on the use of Canvas LMS, Jupyter Notebooks, and a technique called inspections.

## Syllabus Excerpt:

Welcome to BME160 !!

This course is designed to teach programming to scientist/engineers in and around the Life Sciences. This is a hands-on experience where we will move through some of the programming basics and into using Python to build useful tools to help with Life Science Research. We will eventually make use of a few biology-specific programming libraries that the bioinformatics community has been building and making available to us.

We will start off with CodeLab - a tool that we can use to quickly learn programming basics, and we will make use of Jupyter - a tool that lets us write program sections (cells) and to describe those cells in english (markdown cells). We will also reformat that cells into python files that can be run from a low-level environment ( the command-line or "shell" ).

Here is a preview of the main sections of the course:

*   Labs 1 and 2 ( first 2 assignments) are getting experience with the Python "class" environment, how to build objects and how to manipulate "strings" ( strings are lines of text, and we use these a lot in bioinformatics !).
*   Lab 3 is all about iterating and flow control ( making programs work and make decisions), and we will do that by building most of the function that exists in Expasy's ProtParam tool.
*   Lab 4: we work with data files and nucleic acids. In this project, we build a tool to measure properties of coding genes, and can start comparing genomes.
*   Lab 5: here we get to build a gene finder for prokaryotes. This is modeled after NCBI's ORF finder. It is also a real algorithm where we get to use much of Python's power.
*   Lab 6: this one uses multiple objects and multiple sets to solve a real research problem. It is a  custom algorithm -- much like those you will be developing for your own work
*   Final Project: this one .. you and a partner get to pick, working with someone on campus who needs a bit of bioinformatics.

This quarter, we will be making use of **Python 3.10 ( or any version beyond 3.7)**  If you are intending on using your personal computing resources you will no doubt need to either download that version of Python, available for Mac, Linux and Windows (python.org), or better yet, download the very convenient **anaconda** version from **anaconda.org ( it includes Python 3.10 and most everything we will use in this course).**

Learning Objectives
===================

The big one .. we are going to learn Python and its use in bioinformatics. More than that though, we are going to practice translating biological problems into a form that we can answer using a computer and some data. 

With the final project, we also get to communicate our designs in both a written form, as a poster or orally. This will be a real problem that you and your team get to research, analyze, design and implement.

Partners and Teams -- yes, this is a big focus in this course. Basically, no one works by themself in the real world, Inspections -- this is how the real world operates. High-performance teams do all of the best work. We're going to practice this every day.

Just to be clear, we just listed seven of the Program Learning objectives for the BMEB program. The other one is ethical reasoning and it gets highlighted in another of my favorite courses.

Video and Audio Recording[](http://localhost:8888/notebooks/Lecture1.ipynb#Video-and-Audio-Recording)
=====================================================================================================

We will be working in the Active Learning Center, which doesn't have a clear way to record lectures. There may be in-class recordings though.

There may be occasional capture of your voice or video image from time to time in this course. Those recordings are intended for the use of the class as a supplemental learning aid.

Class format
------------

We will be working in a format that maximizes this precious in-person resource that we now have. Class time will include a bit of lecture, a bit of in-class work, and a focus on the weekly-project. To make this work, we have short videos that will provide the essence of our weekly focus. These will need to be viewed before you get to class -- this will let us maximize our in-person time.

We will also be working in small groups throughout the course. I have two intentions here: first, this provides a way for us to learn style and technique from each other; and second, it provides a way for you to debug your programs faster. 

### Inspections

In the commuter industry ( yes, I used to design computer systems ), we use a technique called "Inspections". This is a simple idea where the "author" of a design or program writes down what the intention of their creation is, using text. Remember those Jupyter notebook cells that are just for text ( markdown cells)? Well, that's what these are for. We write down and then present to our group how our stuff works. The group makes suggestions and we then get to save a bunch of time by using the collective wisdom of our team.

### Canvas

All of the assignments and course resources will be available through Canvas. Codelab works through Canvas too ! I will put our datasets on Canvas and our quizzes will be there also.

### Jupyter

We will be starting off using Jupyter. This is a notebook tool that allows embedded python. Please install the latest version of **Anaconda** from **anaconda.com**. Make sure to install the python 3.10 ( or the latest in the Python 3 series) version. Jupyter is built into Anaconda,, along with Numpy, Biopython and lots of other useful modules.

