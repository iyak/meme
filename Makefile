TARGET = meme
all:
	@cd src && $(MAKE)
	@mv src/$(TARGET) .
debug d:
	@cd src && $(MAKE) d
	@mv src/$(TARGET) .
profile prof p:
	@cd src && $(MAKE) p
	@mv src/$(TARGET) .
clean c:
	@cd src && $(MAKE) c
	@rm $(TARGET)

