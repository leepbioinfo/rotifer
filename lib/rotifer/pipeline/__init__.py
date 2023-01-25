class Annotatable:
    def __init__(*args, **kwargs):
        self.results = {}

    def annotate(self,target = 'pid', pipeline='rpsblast', inplace=False, *args, **kwargs):
        '''
        fdsfsd
        '''
        import rotifer.devel.alpha.gian_func as gf
        method = getattr(gf,pipeline)
        if inplace:
            method(self, inplace=inplace, *args, **kwargs)
        else:
            return method(self,*args, **kwargs)
