PROGRAM main
  
  USE mo_message, ONLY: message, message_text
  use mo_ansi_colors, only: color, c_green

  IMPLICIT NONE
  
  Write(*,*) ''
  Write(*,*) 'Test mo_message.f90'

  call message(' mo_message',advance='no')
  message_text = color('o.k.', c_green)
  call message(' ', trim(message_text))

END PROGRAM main
