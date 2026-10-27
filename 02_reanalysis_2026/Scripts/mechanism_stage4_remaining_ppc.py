"""Reuse the validated generic stage-four PPC without the H2/H3-only comparison."""
import inspect
import mechanism_stage4_ppc as base
source=inspect.getsource(base.run)
begin=source.index("        h2=m.ROOT/")
end=source.index("        m.write_json(status,dict(status='complete'",begin)
source=source[:begin]+source[end:]
exec(compile(source,'<stage4_generic_ppc>','exec'),base.__dict__)
if __name__=='__main__':base.run()
