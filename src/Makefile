LIBNAME=cluster
LD=ld

OBJS=	clara.o daisy.o fanny.o meet.o mona.o pam.o spannel.o twins.o

$(LIBNAME).so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so
