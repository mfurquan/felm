module mesh
!========================================================================
! Element, node and zone types
! mesh read/write subroutines
!------------------------------------------------------------------------
! M. Furquan | Feb 2016
!========================================================================
   use global
   implicit none

   type :: node
      real(kind=wp),allocatable,private :: x(nsd)
      integer,private :: id
      contains
         procedure :: set_x
         procedure :: set_id

         procedure :: get_x
         procedure :: get_id
   end type node

   type :: cell
      integer,allocatable,private :: ien(:)
      character(len=celltype_len),private :: celltype
      contains
         procedure :: nen
         procedure :: set_cell

         procedure :: get_celltype
         procedure :: get_ien
   end type cell

   type,extends(cell) :: bface
      real(kind=wp) :: bval(ndf)
      logical :: isDirich=.FALSE.
      contains
         procedure :: set_bval
         procedure :: set_dirich

         procedure :: get_bval
         procedure :: is_dirich
   end type bface

   type :: zone
      integer :: nn,ne,bnf
      type(node),allocatable :: nod(:)
      type(cell),allocatable :: cells(:)
      type(bface),allocatable :: bndry(:)
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

!**************** cell-bound procedures ***********************
     pure function nen(this)
         implicit none
         class(cell),intent(in) :: this
         integer :: nen
         read(this%celltype(celltype_len-1:celltype_len),'(I2)') nen
      end function nen

      subroutine set_cell(this,typname,conn)
         implicit none
         class(cell),intent(inout) :: this
         character(len=celltype_len),intent(in) :: typname
         integer :: conn(:)

         this%celltype = typname
         this%ien = conn
      end subroutine set_cell

      pure function get_celltype(this)
         implicit none
         class(cell),intent(in) :: this
         character(len=celltype_len) :: get_celltype

         get_celltype=this%celltype
      end function get_celltype

      elemental function get_ien(this,inod)
         implicit none
         class(cell),intent(in) :: this
         integer,intent(in) :: inod
         integer :: get_ien

         get_ien = this%ien(inod)
      end function get_ien

!***************** bface-bound procedure ***************************
      subroutine set_bval(this,val)
         class(bface),intent(inout) :: this
         real(kind=wp),intent(in) :: val(:)

         this%bval=val
      end subroutine set_bval

      subroutine set_dirich(this,stat)
         class(bface),intent(inout) :: this
         logical,intent(in) :: stat

         this%isDirich=stat
      end subroutine set_dirich

      elemental function get_bval(this,idf)
         class(bface),intent(in) :: this
         integer,intent(in) :: idf
         real(kind=wp) :: get_bval

         get_bval=this%bval(idf)
      end function get_bval

      elemental function is_dirich(this)
         class(bface),intent(in) :: this
         logical :: is_dirich

         is_dirich=this%isDirich
      end function is_dirich
!***************** zone-bound procedure ****************************
      subroutine read_zone(newzone,nodfile,ienfile,bndfile)
         implicit none
         character(len=fnam_max),intent(in) :: nodfile,ienfile,bndfile
         integer :: Znn,Zne,Zbnf
         class(zone),intent(inout) :: newzone

         real(kind=wp) :: r(nsd),val(ndf)
         integer :: i,j,k,n,ienblock(nen_max)
         character(len=celltype_len) :: tname
         character(len=4) :: fmstr
         logical :: btype

         write(fmstr,'(A2,I1,A)') '(A',celltype_len,')'

         open(10,file=trim(nodfile))
         open(20,file=trim(ienfile))
         read(10,*) Znn
         read(20,*) Zne

         newzone%nn=Znn
         newzone%ne=Zne
         allocate(newzone%nod(Znn))
         allocate(newzone%cells(Zne))

         do i=1,Znn
            read(10,*) (r(j),j=1,nsd)
            call newzone%nod(i)%set_x(r)
         end do
         close(10)

         k=0
         do i=1,Zne
            read(20,fmstr) tname
            read(tname(celltype_len-1:celltype_len),'(I2)') n
            read(20,*) (ienblock(j),j=1,n)
            call newzone%cells(i)%set_cell(tname,ienblock(1:n))
            k=k+n
         end do
         close(10)
         close(20)

         open(30,file=trim(bndfile))
         read(30,*) Zbnf

         newzone%bnf=Zbnf
         allocate(newzone%bndry(Zbnf))

         do i=1,Zbnf
            read(30,fmstr) tname
            read(tname(celltype_len-1:celltype_len),'(I2)') n
            read(30,*) (ienblock(j),j=1,n),btype,(val(j),j=1,ndf)
            call newzone%bndry(i)%set_cell(tname,ienblock(1:n))
            call newzone%bndry(i)%set_bval(val)
            call newzone%bndry(i)%set_dirich(btype)
         end do
      end subroutine read_zone

      subroutine print_zone(pzone,outfile)
         implicit none
         class(zone),intent(in) :: pzone
         character(len=fnam_max),intent(in) :: outfile

         integer :: i,j

         open(10,file=trim(outfile))
         write(10,*) "Elements:"
         do i=1,pzone%ne
            write(10,*) pzone%cells(i)%get_celltype(), &
               (pzone%cells(i)%get_ien(j),j=1,pzone%cells(i)%nen())
         end do
         write(10,*) "Nodes:"
         do i=1,pzone%nn
            write(10,*) (pzone%nod(i)%get_x(j),j=1,nsd)
         end do
         write(10,*) "Boundary faces:"
         do i=1,pzone%bnf
            write(10,*) pzone%bndry(i)%get_celltype(), &
               (pzone%bndry(i)%get_ien(j),j=1,pzone%bndry(i)%nen()), &
               (pzone%bndry(i)%get_bval(j),j=1,ndf),pzone%bndry(i)%is_dirich()
         end do
         close(10)
      end subroutine print_zone
end module mesh
