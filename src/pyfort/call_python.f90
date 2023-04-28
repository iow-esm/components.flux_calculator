module call_python
  use, intrinsic :: iso_c_binding
  implicit none
  
  private

  public :: get
  public :: set
  public :: call_function

  contains

  subroutine check(ret)
    integer :: ret
    if (ret /= 0) stop -1
  end subroutine check

  subroutine set(variable_name, variable_value)
     interface
        function set_(variable_name_c, variable_value_c) result(y) bind (c, name='set')
           use iso_c_binding
           character(c_char) :: variable_name_c 
           integer(c_int) :: variable_value_c
           integer(c_int) :: y
        end function set_
     end interface
     character(len=*) :: variable_name
     integer(c_int) :: variable_value
     character(kind=c_char, len=128) :: variable_name_c
     integer(c_int) :: variable_value_c

     variable_name_c = trim(variable_name)//char(0)
     variable_value_c = variable_value

     call check(set_(variable_name_c, variable_value_c))
  end subroutine set

  subroutine get(variable_name, variable_value)
     interface
        function get_(variable_name_c, variable_value_c) result(y) bind (c, name='get')
           use iso_c_binding
           character(c_char) :: variable_name_c
           integer(c_int) :: variable_value_c(1)
           integer(c_int) :: y
        end function get_
     end interface
     character(len=*) :: variable_name
     integer(c_int) :: variable_value(:)
     character(kind=c_char, len=128) :: variable_name_c
     integer(c_int) :: variable_value_c(size(variable_value, 1))

     variable_name_c = trim(variable_name)//char(0)
     variable_value_c = variable_value

     call check(get_(variable_name_c, variable_value_c))
     variable_value = variable_value_c
  end subroutine get

  subroutine call_function(module_name, function_name)
    interface
       function call_function_(mod_name_c, fun_name_c) result(y) bind(c, name='call_function')
         use iso_c_binding
         character(c_char) mod_name_c, fun_name_c
         integer(c_int) :: y
       end function call_function_
    end interface

    character(len=*) :: module_name, function_name
    character(kind=c_char, len=256) :: mod_name_c, fun_name_c

    mod_name_c = trim(module_name)//char(0)
    fun_name_c = trim(function_name)//char(0)

    call check(call_function_(mod_name_c, fun_name_c))

  end subroutine call_function

end module call_python
