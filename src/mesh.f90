module mesh_TYPE
   use param
   implicit none

   type :: mesh
      integer :: nsd, nnod, nbndry
      real(kind=rp),     allocatable :: x(:,:)
      type(connectivity),allocatable :: bndry(:)
      type(connectivity) :: domain
   contains
      procedure :: read_mesh
   end type mesh

   type :: connectivity
      integer :: ncell, nver
      integer,allocatable :: cell(:,:)
   contains
      procedure :: extract_region
   end type connectivity

   ! nodes per element for Gmsh element types:
   integer,parameter :: tnen(14) =                                       &
      [2,3,4,4,8,6,5,3,6,9,10,27,18,14]

contains
   subroutine read_mesh(this,file_name,domain_id,bndry_id,nsd)
      class(mesh),     intent(in out) :: this
      integer,         intent(in)     :: domain_id, bndry_id(:), nsd
      character(len=*),intent(in)     :: file_name
      integer,allocatable :: elty(:), reg(:), ien(:,:)
      integer             :: file_hande, ne, i, j, n, m, p, eltype
      character           :: dummy

      this%nsd = nsd

      call newunit(file_handle)
      open(file_handle,file=file_name)

      ! read file version
      read(file_handle,*) dummy
      read(file_handle,*) dummy
      if(dummy.llt.'2') error stop "read_mesh: unsupported format"
      read(file_handle,*) dummy

      ! read node coordinates
      read(file_handle,*) dummy
      read(file_handle,*) this%nnod
      do j = 1,this%nnod
         read(file_handle,*) (x(i,j),i = 1,nsd)
      end do
      read(file_handle,*) dummy

      ! read connectivity
      read(file_handle,*) dummy
      read(file_handle,*) ne
      ! temp. arrays for gmsh data
      allocate(elty(ne))
      allocate(reg(ne))
      allocate(ien(nen_max,ne))
      ! read gmsh element section
      do i = 1,ne
         read(file_handle,*) n,elty(i),m,reg(i),p,                      &
            (ien(j,i),j=1,tnen(elty(i)))
      end do
      ! allocate boundary
      this%nbndry = SIZE(bndry)
      allocate(this%bndry(this%nbndry))

      ! initialize domain and boundary meshes
      n = num_vertex(domain_id)
      call domain%extract_region(ien,reg,n,domain_id,ne)
      do i = 1,this%nbndry
         n = num_vertex(bndry_id(i))
         call this%bndry(i)%extract_region(ien,reg,n,bndry_id(i),ne)
      end do

      deallocate(ien)
      deallocate(reg)
      deallocate(elty)
      read(file_handle,*) dummy
   contains
      pure function num_vertex(k)
         integer,intent(in) :: k
         integer :: num_vertex, i
         do i = 1,ne
            if(reg(i)==k) exit
         end do
         num_vertex = tnen(elty(i))
      end function num_vertex
   end subroutine read_mesh

   pure subroutine extract_region(this,ien,reg,nvertex,region_id,ne)
      implicit none
      class(connectivity),intent(in out) :: this
      integer,            intent(in)     :: region_id, ne, nvertex       &
                                            ien(nen_max,ne), reg(ne)
      integer :: m, n

      this%nver  = nvertex
      this%ncell = COUNT(reg == region_id)
      allocate(this%cell(this%nver,this%ncell))

      m = 0
      n = this%nver
      do i = 1,ne; if(reg(i) == region_id) then
         m = m + 1
         this%cell(:,m) = ien(1:n,i)
      end if; end do
   end subroutine extract_region
end module mesh_TYPE
