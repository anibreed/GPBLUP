module link_list
 implicit none
 public :: InsertL, GetL, PrintL, DeleteL

 type,public:: Info
       integer:: id
 end type Info

 type,public:: Node
       type(Info):: Info
       type(Node),pointer:: Next_NP => Null( )
 end type Node

 type(Node),pointer,private:: Root_NP => Null( )

 contains
 recursive subroutine InsertL(Arg_NP, Item)
 type(Node),pointer:: Arg_NP
 type(Info),intent(in):: Item
 if (associated(Arg_NP)) then
    if (Item%id == Arg_NP%Info%id) then
        return
    else if (Arg_NP%Info%id < Item%id) then
        call Insert_Target(Arg_NP,Item)
    else
        call InsertL(Arg_NP%Next_NP,Item)
    end if
 else
    call Insert_Target(Arg_NP,Item )
 end if
 return
 end subroutine

 subroutine Insert_Target(Arg_NP,Item)
 type(Node),pointer:: Arg_NP,Temp_NP
 type(Info),intent(in):: Item
 allocate(Temp_NP)
 Temp_NP%Info = Item
 Temp_NP%Next_NP => Arg_NP
 Arg_NP => Temp_NP
 return
 end subroutine Insert_Target

 function GetL(Arg_NP) result(Item)
 type(Node),pointer:: Arg_NP
 type(Info):: Item
 Item%id=0
 if(associated(Arg_NP)) Item = Arg_NP%Info
 return
 end function GetL

 subroutine DeleteL(Arg_NP)
 type(Node),pointer:: Arg_NP, Temp_NP
 type(Info):: Item
 if (.not.associated(Arg_NP)) return
 Temp_NP => Arg_NP
 Item = Arg_NP%Info
 Arg_NP => Arg_NP%Next_NP
 deallocate(Temp_NP)
 end subroutine

 recursive subroutine PrintL(Arg_NP)
 type (Node),pointer:: Arg_NP, Trav_NP
 if (associated(Arg_NP)) then
     call Print_Target(Arg_NP)
     call PrintL(Arg_NP % Next_NP)  ! Advance to next node
 end if
 return
 end subroutine PrintL

 subroutine Print_Target(Arg_NP)
 type(Node),pointer:: Arg_NP
 print *, " Print: ", Arg_NP%Info
 return
 end subroutine Print_Target

end module link_list

