#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct
{
    unsigned char b;
    unsigned char g;
    unsigned char r;
}pixel;
typedef struct
{
    unsigned char* header;
    pixel* pixels;
    unsigned int width;
    unsigned int height;
}imagine;
imagine get_img(char* path)
{
    imagine A;
    FILE* f = fopen(path , "rb");
    A.header = malloc( 54*sizeof(char) );
    fread(A.header , sizeof(char) , 54 , f);
    A.width =  *(unsigned int *)&A.header[18];
    A.height =  *(unsigned int *)&A.header[22];
    A.pixels = malloc(A.width*A.height*sizeof(pixel));
    unsigned int octpad = (4 - ((3*A.width) % 4)) % 4 ;
    unsigned char a;
    for(int i=A.height-1 ; i>=0 ; i--)
    {
        for (unsigned int j = 0; j < A.width; j++)
        {
            fread(&A.pixels[i * A.width + j].b, sizeof(char), 1, f);
            fread(&A.pixels[i * A.width + j].g, sizeof(char), 1, f);
            fread(&A.pixels[i * A.width + j].r, sizeof(char), 1, f);
        }
        if(octpad>0)
            for(int j = 0 ; j < octpad ; j++)
                fread(&a, sizeof(char), 1, f);
    }
    fclose(f);
    return A;
}
void put_img(imagine A , char* path_out)
{
    FILE* g = fopen(path_out , "w+b");
    fwrite(A.header , sizeof(char) , 54 , g);
    unsigned int octpad = (4 - ((3*A.width) % 4)) % 4 ;
    unsigned char a = 0;
    for(int i=A.height-1 ; i>=0 ; i--)
    {
        for (unsigned int j = 0; j < A.width; j++)
        {
            fwrite(&A.pixels[i * A.width + j].b, 1, 1, g);
            fwrite(&A.pixels[i * A.width + j].g, 1, 1, g);
            fwrite(&A.pixels[i * A.width + j].r, 1, 1, g);
        }
        if(octpad>0)
            for(int j=0 ; j<octpad ; j++)
                fwrite(&a, 1, 1, g);
    }
    fclose(g);
}

//********************************************MODUL CRIPTARE********************************************

unsigned int xorshift32(unsigned int seed)
{
    unsigned int x = seed;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    seed = x;
    return seed;
}
pixel pixels_xor(pixel P, pixel Q)
{
    P.b = P.b ^ Q.b;
    P.g = P.g ^ Q.g;
    P.r = P.r ^ Q.r;
    return P;
}
pixel pixel_xor_int(pixel P, unsigned int x)
{
    unsigned int aux1 = x % 256;
    P.b = P.b ^ (unsigned char)aux1;
    x = x / 256;
    unsigned int aux2 = x % 256;
    P.g = P.g ^ (unsigned char)aux2;
    x = x / 256;
    unsigned int aux3 = x % 256;
    P.r = P.r ^ (unsigned char)aux3;
    return P;
}
unsigned int* generate_sequence(unsigned int seed, unsigned int width, unsigned int height)
{
    unsigned int* v = malloc(2*width*height*sizeof(int));
    v[0] = seed;
    for(unsigned int i=1 ; i<2*width*height ; i++)
        v[i] = xorshift32(v[i-1]);
    return v;
}
unsigned int* Durstenfeld(unsigned int seed, unsigned int width, unsigned int height)
{
    unsigned int* p = malloc(width*height*sizeof(int));
    for(unsigned int i=0 ; i<width*height ; i++)
        p[i]=i;
    unsigned int* v;
    v = generate_sequence(seed, width, height);
    int j=1, k = width*height;
    for(int i=width*height-1; i>=1; i--)
    {
        unsigned int aux = p[v[j] % k];
        p[v[j] % k] = p[i];
        p[i] = aux;
        j++;
        k--;
    }
    free(v);
    return p;
}
void mix_pixels(imagine A , unsigned int seed)
{
    unsigned int* p = Durstenfeld(seed , A.width , A.height);
    put_img(A , "auxiliar.bmp");
    imagine B = get_img("auxiliar.bmp");
    for(unsigned int i=0 ; i<A.width*A.height ; i++)
        A.pixels[p[i]] = B.pixels[i];
    free(p);
}
void crypt_image(char* path_initial, char* path_final, char* path_key)
{
    unsigned int seed , SV;
    FILE* f = fopen(path_key , "r");
    fscanf(f,"%u %u",&seed,&SV);
    fclose(f);
    imagine A = get_img(path_initial);
    unsigned int* v = generate_sequence(seed , A.width , A.height);
    mix_pixels(A , seed);
    A.pixels[0] = pixel_xor_int(pixel_xor_int(A.pixels[0] , SV) , v[A.width*A.height]);
    for(unsigned int i=1 ; i<A.width*A.height ; i++)
        A.pixels[i] = pixel_xor_int(pixels_xor(A.pixels[i-1] , A.pixels[i]) , v[A.width*A.height+i]);
    put_img(A , path_final);
    free(v);
}
unsigned int* reverse_Durstenfeld(unsigned int seed, unsigned int width, unsigned int height)
{
    unsigned int* v = Durstenfeld(seed , width , height);
    unsigned int* w = malloc(width*height*sizeof(int));
    unsigned int aux;
    for(unsigned int i=0 ; i<width*height ; i++)
        w[v[i]] = i;
    free(v);
    return w;
}
void unmix_pixels(imagine A, unsigned int seed)
{
    unsigned int* r = reverse_Durstenfeld(seed , A.width , A.height);
    put_img(A , "auxiliar.bmp");
    imagine B = get_img("auxiliar.bmp");
    for(unsigned int i=0 ; i<A.width*A.height ; i++)
        A.pixels[r[i]] = B.pixels[i];
    free(r);
}
void decrypt_image(char* path_initial, char* path_final, char* path_key)
{
    unsigned int seed , SV;
    FILE* f = fopen(path_key , "r");
    fscanf(f,"%u %u",&seed,&SV);
    fclose(f);
    imagine A = get_img(path_initial);
    unsigned int* v = generate_sequence(seed , A.width , A.height);
    for(int i=A.width*A.height-1 ; i>0 ; i--)
        A.pixels[i] = pixel_xor_int(pixels_xor(A.pixels[i] , A.pixels[i-1]) , v[A.width*A.height+i]);
    A.pixels[0] = pixel_xor_int(pixel_xor_int(A.pixels[0] , SV) , v[A.width*A.height]);
    unmix_pixels(A , seed);
    put_img(A , path_final);
    free(v);
}
typedef struct
{
    double b;
    double g;
    double r;
}chi_values;
chi_values chi_squared(char* path)
{
    imagine A = get_img(path);
    chi_values C;
    double f = (double)A.width*(double)A.height/256;
    //blue chi
    double* f_b = malloc(256*sizeof(double));
    for(unsigned int i=0 ; i<256 ; i++)  f_b[i]=0;
    for(unsigned int i=0 ; i<256 ; i++)
        for(unsigned int j=0 ; j<A.height*A.width ; j++)
            if(i == A.pixels[j].b) f_b[i]++;
    C.b = 0;
    for(unsigned int i=0 ; i<256 ; i++)
        C.b += (f_b[i] - f)*(f_b[i] - f) / f;
    //green chi
    double* f_g = malloc(256*sizeof(double));
    for(unsigned int i=0 ; i<256 ; i++)  f_g[i]=0;
    for(unsigned int i=0 ; i<256 ; i++)
        for(unsigned int j=0 ; j<A.height*A.width ; j++)
            if(i == A.pixels[j].g) f_g[i]++;
    C.g = 0;
    for(unsigned int i=0 ; i<256 ; i++)
        C.g += (f_g[i] - f)*(f_g[i] - f) / f;
    //red chi
    double* f_r = malloc(256*sizeof(double));
    for(unsigned int i=0 ; i<256 ; i++)  f_r[i]=0;
    for(unsigned int i=0 ; i<256 ; i++)
        for(unsigned int j=0 ; j<A.height*A.width ; j++)
            if(i == A.pixels[j].r) f_r[i]++;
    C.r = 0;
    for(unsigned int i=0 ; i<256 ; i++)
        C.r += (f_r[i] - f)*(f_r[i] - f) / f;

    return C;
}

//***************************************MODUL RECUNOASTERE PATTERNURI***************************************
void grayscale(imagine A)
{
    unsigned char gray;
    for(unsigned int i=0 ; i<=A.height*A.width ; i++)
    {
        gray = (unsigned char)(0.299 * A.pixels[i].r + 0.587 * A.pixels[i].g + 0.114 * A.pixels[i].b);
        A.pixels[i].r = gray;
        A.pixels[i].g = gray;
        A.pixels[i].b = gray;
    }
}
typedef struct
{
    unsigned int corner;
    double match;
    pixel color;
}detections;
void template_match(detections**T , imagine A , imagine S , pixel C , float p , int* matches)
{
    double corr = 0;
    int gasite_aici = 0;
    for(unsigned int i=0 ; i<A.height * A.width ; i++)
    {
        if ((i / A.width < A.height - S.height) && (i % A.width < A.width - S.width)) {
            double f = 0, s = 0, sigmaf = 0, sigmas = 0;
            for (unsigned int j = 0; j < S.height * S.width; j++) {
                f = f + A.pixels[i + (j / S.width) * A.width + j % S.width].b;
                s = s + S.pixels[j].b;
            }
            f = f / (S.width * S.height);
            s = s / (S.width * S.height);
            for (unsigned int j = 0; j < S.height * S.width; j++) {
                sigmaf = sigmaf + (A.pixels[i + (j / S.width) * A.width + j % S.width].b - f) *
                                  (A.pixels[i + (j / S.width) * A.width + j % S.width].b - f) /
                                  (S.width * S.height - 1);
                sigmas = sigmas + (S.pixels[j].b - s) * (S.pixels[j].b - s) / (S.width * S.height - 1);
            }
            sigmaf = sqrt(sigmaf);
            sigmas = sqrt(sigmas);
            for (unsigned int j = 0; j < S.height * S.width; j++) {
                corr = corr + (A.pixels[i + (j / S.width) * A.width + j % S.width].b - f) * (S.pixels[j].b - s);
            }
            corr = corr / (S.height * S.width);
            corr = corr / (sigmaf * sigmas);
            if (corr > p && corr < 1) {
                (*T) = realloc((*T), ((*matches) + 1) * sizeof(detections));
                (*T)[(*matches)].corner = i;
                (*T)[(*matches)].match = corr;
                (*T)[(*matches)].color = C;
                (*matches)++;
                gasite_aici++;
            }
            corr = 0;
        }
    }

    //printf("Gasite aici: %d\n", gasite_aici);
}
void draw_window(imagine A , detections t , unsigned int width , unsigned int height)
{
    for(unsigned int i = 0; i < width ; i++)
    {
        A.pixels[t.corner + i] = t.color;
        A.pixels[t.corner + A.width* height + i] = t.color;
    }
    for(int i = 0; i < height ; i++)
    {
        A.pixels[t.corner + i*A.width] = t.color;
        A.pixels[t.corner + width + i*A.width] = t.color;
    }
}
int cmp(const void *a , const void *b)
{
    if(((detections *)a)->match > ((detections *)b)->match)
        return  1;
    return -1;
}
unsigned int area_intersection(imagine A , unsigned int a  , unsigned int b , unsigned int H , unsigned int W)
{
    int Area,area;
    Area = W*H*2;
    int h,w;
    h = a/A.width - b/A.width;
    if(h < 0) h = -h;
    w = a % A.width - b % A.width;
    if(w < 0) w = -w;
    if(w<W && h<H && (H-h)*(W-w)>0.2*Area)
        return 1;
    return 0;
}
void remove_non_maximum(detections* T , imagine A , imagine S , int N)
{
    for(int i=0; i<N-1 ; i++)
    {
        if(T[i].match)
            for(int j=i ; j<N ; j++)
            {
                if(area_intersection(A , T[i].corner , T[j].corner , S.height , S.width))
                    T[j].match = 0;
            }
    }
}

int main()
{
    //CRIPTARE**************//
    char initial[50] , crypt[50] , key[50] , crypted[50] , decrypted[50] , key2[50] ,da=2;
    while(da!=0&&da!=1)
    {
        printf("Pt criptare dati 1; pt fara criptare dati 0")
        scanf("%u",&da);
    }
    if(da==1)
    {
        printf("Calea imaginii initiale este ");
        scanf("%s", initial);
        printf("\nCalea la care se face criptarea este ");
        scanf("%s", crypt);
        printf("\nCalea cheii este ");
        scanf("%s", key);
        crypt_image(initial, crypt, key);

        printf("Calea imaginii criptate este ");
        scanf("%s", crypted);
        printf("\nCalea la care se face decriptarea este ");
        scanf("%s", decrypted);
        printf("\nCalea cheii este ");
        scanf("%s", key2);
        decrypt_image(crypted, decrypted, key2);


        chi_values original_chi = chi_squared(initial);
        printf("Valorile testului chi pentru imaginea initiala pe canalele de culoare R, G, respectiv B sunt %f %f %f\n", original_chi.r, original_chi.g, original_chi.b);
        chi_values crypted_chi = chi_squared(crypt);
        printf("Valorile testului chi pentru imaginea criptata pe canalele de culoare R, G, respectiv B sunt %f %f %f", crypted_chi.r, crypted_chi.g, crypted_chi.b);
    }

    //PATTERNURI***********//
    printf("Calea imaginii initiale este ");
    scanf("%s",initial);
    imagine A = get_img(initial);
    grayscale(A);
    imagine S;

    detections* D = malloc(sizeof(detections));
    int matches = 0 , sab;

    printf("Cate sabloane bagam?");
    scanf("%u",&sab);
    FILE *f = fopen("pixel_colors.txt","r");
    for (int j = 0; j < sab ; ++j)
    {
        printf("Sablon %d:", j + 1);
        scanf("%s",key);
        S = get_img(key);
        grayscale(S);
        pixel C;
        fscanf(f, "%u", &C.b); fscanf(f, "%u", &C.g); fscanf(f, "%u", &C.r);

        template_match(&D , A , S , C , 0.5, &matches);
    }

    qsort(D , matches , sizeof(detections) , cmp);
    remove_non_maximum(D,A,S,matches);
    qsort(D , matches , sizeof(detections) , cmp);

    for (int i = 0; i < matches; ++i)
    {
        if (D[i].match)
        {
            draw_window(A, D[i] , S.width , S.height);
        }

    }
    put_img(A,"result.bmp");

    return 0;
}