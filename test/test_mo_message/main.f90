PROGRAM main
  
  USE mo_message, ONLY: message, message_text

  IMPLICIT NONE
  
  Write(*,*) ''
  Write(*,*) 'Test mo_message.f90'

  call message(' mo_message',advance='no')
  message_text = 'o.k.'
  call message(' ', trim(message_text))

END PROGRAM main
