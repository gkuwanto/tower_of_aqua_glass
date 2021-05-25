import pygame
import numpy as np
from map import Map

pygame.init()
WIDTH = 800
HEIGHT = 600
SCALE = 10

gameDisplay = pygame.display.set_mode((WIDTH,HEIGHT))

pygame.display.set_caption('Menara Gelas')
clock = pygame.time.Clock()

def to_game_coord(math_coord):
    x,y = math_coord
    return (SCALE* x), (HEIGHT - (SCALE* y) - 30)
def to_math_coord(game_coord):
    x,y = game_coord
    return x / SCALE, (HEIGHT-y)/SCALE



# Gelas Local
def get_pentagon(x,y):
    return [(x-22.5, y-97), (x+22.5,y-97), (x+37.5, y-2), (x,y), (x-37.5,y-2)]
def get_level(mouse_x, gelas):
    cnt = 0
    for x,y, below in gelas:
        left = x - 75
        right = x + 75
        if left < mouse_x < right:
            cnt = max(cnt, below+1)
    return cnt

m = Map()
closed = False
reset = False
gelas = []
level_to_upper_interval = {-1: []}
last_stable = 2

while not closed:
    if reset == True:
        del m
        del gelas
        del level_to_upper_interval
        m = Map()
        reset = False
        gelas = []
        level_to_upper_interval = {-1: []}
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            closed = True
        if event.type == pygame.MOUSEMOTION:
            mouse_x, mouse_y = pygame.mouse.get_pos()
            below = get_level(mouse_x, gelas)
            mouse_y = 590 - (below * 97)
        if event.type == pygame.MOUSEBUTTONDOWN:
            gelas.append((mouse_x,mouse_y, below))
            try:
                level_to_upper_interval[below].append((mouse_x-22.5, mouse_x+22.5))
            except:
                level_to_upper_interval[below] = [(mouse_x-22.5, mouse_x + 22.5)]
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_RETURN:
                if len(gelas):
                    level_to_index = {}
                    idx_p = 0
                    for gx, gy, level in sorted(gelas, key=lambda x: x[2]):
                        x, y = to_math_coord((gx,gy))
                        try:
                            level_to_index[level].append(idx_p)
                        except:
                            level_to_index[level] = [idx_p]
                        tahan_kiri = False
                        tahan_kanan = False
                        for min_bound_x, max_bound_x in level_to_upper_interval[level-1]:
                            if min_bound_x< gx-37.5 <max_bound_x:
                                tahan_kiri = True
                            if min_bound_x< gx+37.5 <max_bound_x:
                                tahan_kanan = True
                            
                        m.add_gelas(x, y, level==0, tahan_kiri, tahan_kanan)
                        for i in range(level-1, 0, -1):
                            for idx in level_to_index[i]:
                                m.add_interaction(idx_p, idx)
                        idx_p += 1
                    m.interact_all()
                    last_stable = m.is_stable()
                    print(m.get_total_center_of_mass())
                    m = Map()
                    print(last_stable)
            if event.key == pygame.K_r:
                reset = True
                last_stable = 2

    
    gameDisplay.fill((255,255,255))
    pygame.draw.polygon(gameDisplay, (200,200,255), get_pentagon(mouse_x,mouse_y))
    for gx, gy, _ in gelas:
        pygame.draw.polygon(gameDisplay, (0,0,255), get_pentagon(gx,gy))
    if last_stable == 2:
        pass
    elif last_stable == 1:
        pygame.draw.circle(gameDisplay,(0,255,0), (60,60),50 )
    else:
        pygame.draw.circle(gameDisplay,(255,0,0), (60,60),50 )

    pygame.display.update()
    clock.tick(60)


pygame.quit()
quit()