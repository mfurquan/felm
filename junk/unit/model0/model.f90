module model
   use global
   implicit none

   type :: node
      real(kind=wp),private :: x(nsd),u(ndf)
      integer,private :: id
      contains
         procedure :: set_x
         procedure :: set_u
         procedure :: set_id

         procedure :: get_x
         procedure :: get_u
         procedure :: get_id
   end type node

   integer,protected,allocatable,target :: ienblock(:)

   type :: element
      integer,pointer,private :: ien(:) => null()
      character(len=elty_len),private :: typ
      contains
         procedure :: nen
         procedure :: set_element

         procedure :: get_typ
         procedure :: get_ien
   end type element

   type :: zone
      integer :: nn,ne
      type(node),allocatable :: nod(:)
      type(element),allocatable :: elem(:)
      contains
         procedure :: read_zone
         procedure :: print_zone
   end type zone

   contains

!********************* node-bound procedures ************************
      subroutine set_x(this,coor)
         implicit none
         class(node),intent(inout) :: this
         real(kind=wp),intent(in) :: coor(nsd)

         this%x=coor
      end subroutine set_x

      subroutine set_u(this,disp)
         implicit none
         class(node),intent(inout) :: this
         real(kind=wp),intent(in) :: disp(ndf)

         this%u=disp
      end subroutine set_u

      subroutine set_id(this,idx)
         implicit none
         class(node),intent(inout) :: this
         integer,intent(in) :: idx

         this%id=idx
      end subroutine set_id

      elemental function get_x(this,icor)
         implicit none
         class(node),intent(in) :: this
         real(kind=wp) :: get_x
         integer,intent(in) :: icor

         get_x=this%x(icor)
      end function get_x

      elemental function get_u(this,idof)
         implicit none
         class(node),intent(in) :: this
         real(kind=wp) :: get_u
         integer,intent(in) :: idof

         get_u=this%u(idof)
      end function get_u

      pure function get_id(this)
         implicit none
         class(node),intent(in) :: this
         integer :: get_id

         get_id=this%id
      end function get_id

!**************** element-bound procedures ***********************
     pure function nen(this)
         implicit none
         class(element),intent(in) :: this
         integer :: nen
         read(this%typ(elty_len-1:elty_len),'(I2)') nen
      end function nen

      subroutine set_element(this,typname,iblock)
         implicit none
         class(element),intent(inout) :: this
         character(len=elty_len) :: typname
         integer :: iblock,nblock

         this%typ=typname
         read(typname(elty_len-1:elty_len),'(I2)') nblock
         this%ien => ienblock(iblock:iblock+nblock-1)
      end subroutine set_element

      pure function get_typ(this)
         implicit none
         class(element),intent(in) :: this
         character(len=elty_len) :: get_typ

         get_typ=this%typ
      end function get_typ

      elemental function get_ien(this,inod)
         implicit none
         class(element),intent(in) :: this
         integer,intent(in) :: inod
         integer :: get_ien

         get_ien = this%ien(inod)
      end function get_ien

!***************** zone-bound procedure ****************************
      subroutine read_zone(newzone,nodfile,ienfile,Znn,Zne,Znne)
         implicit none
         character(len=fnam_max),intent(in) :: nodfile,ienfile
         integer,intent(in) :: Znn,Zne,Znne
         class(zone),intent(inout) :: newzone

         real(kind=wp) :: r(nsd)
         integer :: i,j,k,n
         character(len=elty_len) :: tname
         character(len=4) :: fmstr

         write(fmstr,'(A2,I1,A)') '(A',elty_len,')'

         newzone%nn=Znn
         newzone%ne=Zne
         allocate(newzone%nod(Znn))
         allocate(newzone%elem(Zne))

         open(10,file=trim(nodfile))
         do i=1,Znn
            read(10,*) (r(j),j=1,nsd)
            call newzone%nod(i)%set_x(r)
         end do
         close(10)

         allocate(ienblock(Znne))
         open(20,file=trim(ienfile))
         k=0
         do i=1,Zne
            read(20,fmstr,advance='no') tname
            read(tname(elty_len-1:elty_len),'(I2)') n
            read(20,*) (ienblock(j),j=k+1,k+n)
            call newzone%elem(i)%set_element(tname,k+1)
            k=k+n
         end do

         if(debugOn) then
            if(k/=Znne) write(*,*) 'Error: Insufficient allocation in ienblock'
            stop
         end if
      end subroutine read_zone

      subroutine print_zone(pzone,outfile)
         implicit none
         class(zone),intent(in) :: pzone
         character(len=fnam_max),intent(in) :: outfile

         integer :: i,j

         open(10,file=trim(outfile))
         do i=1,pzone%ne
            write(10,*) (pzone%elem(i)%get_ien(j),j=1,pzone%elem(i)%nen())
         end do
         do i=1,pzone%nn
            write(10,*) (pzone%nod(i)%get_x(j),j=1,nsd)
         end do
         close(10)
      end subroutine print_zone
end module model
