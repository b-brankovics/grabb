all:		compile moveExec
	
compile:	
		@cd ./src; make

moveExec:	
		@cp ./src/edena ./bin/.

clean:		
		@cd ./src; make clean
		@cd ./src; rm edena
		@cd ./bin; rm edena
	
