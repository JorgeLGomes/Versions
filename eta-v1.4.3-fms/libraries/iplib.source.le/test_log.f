        program test

        logical*1 li(25)

        do i=1,25

        if (mod(i,2) .eq. 0) then
        li(i)=.true.
        else
        li(i)=.false.
        endif

        enddo

        write(6,*) 'li= ', (li(i),i=1,25)
        end
