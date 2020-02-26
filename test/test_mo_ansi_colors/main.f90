program test_ansi_colors

  use mo_ansi_colors, only: color, c_red, c_green, c_yellow, c_blue, c_magenta, c_cyan, c_black, c_white

  implicit none

  ! character(len=*), parameter :: endl = new_line('a')

  ! print '(a)', &
  !      color('Red',     c_red)     // endl // &
  !      color('Green',   c_green)   // endl // &
  !      color('Yellow',  c_yellow)  // endl // &
  !      color('Blue',    c_blue)    // endl // &
  !      color('Magenta', c_magenta) // endl // &
  !      color('Cyan',    c_cyan)    // endl // &
  !      color('Black',   c_black)   // endl // &
  !      color('White',   c_white)

  write(*,*) 'This is ', color('Red',     c_red),     '.'
  write(*,*) 'This is ', color('Green',   c_green),   '.'
  write(*,*) 'This is ', color('Yellow',  c_yellow),  '.'
  write(*,*) 'This is ', color('Blue',    c_blue),    '.'
  write(*,*) 'This is ', color('Magenta', c_magenta), '.'
  write(*,*) 'This is ', color('Cyan',    c_cyan),    '.'
  write(*,*) 'This is ', color('Black',   c_black),   '.'
  write(*,*) 'This is ', color('White',   c_white),   '.'
  
  write(*,*) ''
  write(*,*) 'mo_ansi_colors ', color('o.k.', c_green)
  
end program test_ansi_colors
