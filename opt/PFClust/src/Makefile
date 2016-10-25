# Makefile for metsite
JAVAC=javac
JAR=jar

all: compileMain buildJar

clean:
	rm clustering/*.class main/*.class utilities/*.class

install:
	mv PFClust.jar ../

compileMain:
	$(JAVAC) main/Main.java

buildJar:
	$(JAR) cvfm PFClust.jar META-INF/MANIFEST.MF clustering/*.class main/*.class utilities/*.class
