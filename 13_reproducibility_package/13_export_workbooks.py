"""Portable XLSX export; cell values come unchanged from the pipeline JSON."""
from pathlib import Path
import json, sys, math, textwrap
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from PIL import Image, ImageDraw, ImageFont
import matplotlib

def main():
    build = Path(sys.argv[1])
    specs = json.loads((build/'workbooks.json').read_text(encoding='utf-8'))
    qa = build/'review_materials/workbook_visual_check'
    qa.mkdir(parents=True, exist_ok=True)
    fontdir = Path(matplotlib.get_data_path())/'fonts/ttf'
    regular = ImageFont.truetype(str(fontdir/'DejaVuSans.ttf'), 14)
    bold = ImageFont.truetype(str(fontdir/'DejaVuSans-Bold.ttf'), 14)
    for rel, sheets in specs.items():
        wb = Workbook(); wb.remove(wb.active)
        for spec in sheets:
            ws = wb.create_sheet(spec['name'])
            for row in spec['rows']: ws.append(row)
            ws.freeze_panes = 'A2'; ws.auto_filter.ref = ws.dimensions
            ws.sheet_view.showGridLines = False
            widths = []
            for i, col in enumerate(ws.iter_cols(), 1):
                width = min(46, max(12, max(len(str(c.value or '')) for c in col)+2))
                widths.append(width)
                ws.column_dimensions[get_column_letter(i)].width = width
                for c in col:
                    c.font = Font(name='DejaVu Sans', size=10)
                    c.alignment = Alignment(vertical='top',wrap_text=True)
                    if isinstance(c.value,float): c.number_format='0.000000'
            for row in ws.iter_rows():
                lines=max(max(1,math.ceil(len(str(c.value or ''))/(widths[i]*0.75))) for i,c in enumerate(row))
                ws.row_dimensions[row[0].row].height=max(28,15*lines+10)
            for c in ws[1]:
                c.fill=PatternFill('solid',fgColor='17365D');c.font=Font(name='DejaVu Sans',bold=True,color='FFFFFF',size=10)
            # Portable previews use the same rows and bundled DejaVu font; no Office renderer.
            rows=spec['rows'][:16]; cols=min(9,len(rows[0])); sizes=[int(w*8) for w in widths[:cols]]
            wrapped=[]; heights=[]
            for row in rows:
                cells=[]
                font=bold if not wrapped else regular
                for i,v in enumerate(row[:cols]):
                    lines=[];line=''
                    for ch in str(v if v is not None else ''):
                        if font.getlength(line+ch)>sizes[i]-12:
                            lines.append(line);line=''
                        line+=ch
                    lines.append(line);cells.append(lines)
                wrapped.append(cells);heights.append(max(len(c) for c in cells)*19+14)
            img=Image.new('RGB',(sum(sizes)+1,sum(heights)+1),'white');draw=ImageDraw.Draw(img);y=0
            for ri,cells in enumerate(wrapped):
                x=0
                for ci,lines in enumerate(cells):
                    draw.rectangle((x,y,x+sizes[ci],y+heights[ri]),fill='#17365D' if ri==0 else ('#F0F4F8' if ri%2==0 else 'white'),outline='#D8DFE6')
                    draw.multiline_text((x+5,y+6),'\n'.join(lines),font=bold if ri==0 else regular,fill='white' if ri==0 else '#172B4D',spacing=3)
                    x+=sizes[ci]
                y+=heights[ri]
            img.save(qa/(Path(rel).stem+'_'+spec['name']+'.png'))
        path=build/rel;path.parent.mkdir(parents=True,exist_ok=True);wb.save(path)
    (build/'review_materials/workbook_data_snapshot.json').write_text(json.dumps(specs,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'OPENPYXL_EXPORT_PASS: {len(specs)} workbooks; {sum(map(len,specs.values()))} sheets')

if __name__=='__main__': main()
