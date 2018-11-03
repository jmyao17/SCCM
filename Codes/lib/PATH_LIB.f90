  character(500) function find_file(directory_name,file_name)
  ! returns directory name from found from environment variable
  character(*) :: file_name,directory_name
  character(1000) :: path
  integer :: istart,iend,i
  logical :: isthere,last_chance

  ! file_name = file_name)
  ! directory_name = adjustl(directory_name)

  call GETENV(trim(directory_name),path)
  path=adjustl(path)

  istart = 1
  do while (.true.)
     i = istart
     do while (path(i:i).ne.':')
        iend = i
        i = i+ 1
        if (path(i:i) == ' ') then
           i = 0
           exit
        end if
     end do
     inquire(file=path(istart:iend)//'/'//trim(file_name),exist=isthere)
     if (isthere) then
        find_file = path(istart:iend)//'/'
        return
     else if (i == 0) then
        inquire(file='./'//file_name,exist=last_chance)
        if (last_chance) then
           find_file = './'
           return
        else
           PRINT*, 'FILE NOT FOUND: ',trim(adjustl(file_name))
           print*, 'CHECK ENVIRONMENT VARIABLE: ',trim(adjustl(directory_name))
           STOP
        end if
     end if
     istart = iend+2
  end do
end function find_file
